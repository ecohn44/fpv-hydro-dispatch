using JuMP
using Gurobi
using MAT
using Printf
using CSV
using DataFrames
using Dates
using Plots
using Base.Filesystem
include("plots.jl")

# ----------------- UNIT CONVERSION ----------------- #

function cfs_to_m3s(cfs)
    # Conversion factor: 1 cfs is approximately 0.0283168 m³/s
    conversion_factor = 0.0283168
    return cfs * conversion_factor
end

function af_to_m3(acft)
    # Conversion factor: 1 acre-foot is approximately 1233.48 cubic meters
    conversion_factor = 1233.48
    return acft * conversion_factor
end

function dailyflow_to_hourly(q, T)
    # Input: Vector q of average daily flow values
    # Output: Hourly flow vector held constant over the day
    hourly_q = []
    for i in q
        rpt = fill(i, T)
        append!(hourly_q, rpt)
    end
    return hourly_q
end

# ----------------- DATA LOAD ----------------- #
# 1 WEEK SIMULATION FOR JAN 1-7 2022
year = "2022";
month = "January";
month_num = 1;

LMP_path = string("data/LMP-meads-2-N101-",month,year,".csv");
netflow_path = string("data/netflow-daily-",year,"-USBR.csv");
solarrad_path = "data/lake-mead-monthly-solar-rad-capacity.csv";
release_path = string("data/mead-totalrelease-",month,year,".csv");
storage_path = string("data/USBR-lakemead-daily-",month,year,".csv");

# Read in Local Marginal Price for Lake Mead
RTP = DataFrame(CSV.File(LMP_path));
RTP_2d = [row.LMP for row in eachrow(RTP)];

# Read in system inflow (lake powell inflow)
netflow = DataFrame(CSV.File(netflow_path));
# convert to m3/s from cfs
netflow_2d = [cfs_to_m3s(row["powell_inflow (cfs)"]) for row in eachrow(netflow)];

# Read in solar radiation capacity
alpha = DataFrame(CSV.File(solarrad_path));
alpha_2d = [row[month] for row in eachrow(alpha)];
min_val = minimum(alpha_2d);
max_val = maximum(alpha_2d);
alpha_norm = (alpha_2d .- min_val) / (max_val - min_val); # Normalize data 

# Read in Water Contract Data
release = DataFrame(CSV.File(release_path));
# convert from ac-ft to m3
release_2d = [af_to_m3(row.Result) for row in eachrow(release)];

# Read in historical storage levels 
storage = DataFrame(CSV.File(storage_path));
# convert from ac-ft to m3
storage_2d = [af_to_m3(row["storage (ac-ft)"]) for row in eachrow(storage)];

# ----------------- SIMULATION PARAMETERS ----------------- #
## Overview: Power contractors operate on weekly subcontracts. 
# Therefore, will let the simulation run every hour of the week before resetting for the next water contract 

T = 24; # time step per day
N_sim = 7; # number of days (per week)
W_sim = 1; # number of weeks 

# Extract data for time period of interest
L = RTP_2d[1:T*N_sim,1];    # hourly LMP [$/MWh]
Q = netflow_2d[1:N_sim,1];  # daily system net flow into [m3]
U = release_2d[1:N_sim,1]   # daily total release contract [m3]
S = storage_2d[1:N_sim,1]   # historical daily storage levels [m3]
q = dailyflow_to_hourly(Q, T); # hourly system net flow (constant over each day) [m3]
alpha_norm_w = repeat(alpha_norm, N_sim)

e = 0    # evaporation constant [m/day]
min_ut = cfs_to_m3s(5000)    # min daily release limit [m3/s]
max_ut = cfs_to_m3s(25000)   # max daily release limit [m3/s]
RR_dn = cfs_to_m3s(-2500) # down ramp rate limit [m3/s]
RR_up = cfs_to_m3s(4000)  # up ramp rate limit [m3/s]
PF = 3300   # max feeder capacity (3.3 GW) 
PS = 1000   # max solar field capacity [MW] (1 GW) 
SA = 1.3e9   # surface area of both reservoirs [m^2] (504 sq miles)
eta = .9     # efficiency of release-energy conversion
rho_w = 1000 # density of water [kg/m^3]
g = 9.8      # acceleration due to gravity [m/s^2]
s2hr = 3600  # seconds in an hour 
ws2MWh = 1e6 # conversion from Ws to MWh

# ----------------- VARIABLES ----------------- #

# Weekly parameters 
Uw = sum(U[1:N_sim]); # Weekly Water contract
S0 = 6.3e10 #S[1]   # initial reservoir conditions [m3]

# initialize optimization model
model = Model(Gurobi.Optimizer)
set_silent(model) # no outputs

@variable(model, ph[1:T*N_sim], lower_bound = 0)      # total hydro discharge power [Mw]
@variable(model, ps[1:T*N_sim], lower_bound = 0)      # total solar discharge power [MW]
@variable(model, S[1:T*N_sim], lower_bound = 0)       # hourly storage level [m3]
@variable(model, u[1:T*N_sim], lower_bound = 0)       # water release [m3/s]
@variable(model, R)                                   # market revenue [$]

# ----------------- OPTIMIZATION  ----------------- #

# operations revenue
@constraint(model, Revenue, R == sum(L.*(ph + ps)))

# (1a) initial reservoir mass balance equation
@constraint(model, MassBalInit, S[1] == S0 + s2hr*(q[1] - u[1])) # .- e*SA)
# (1b) rest reservoir mass balance equation
@constraint(model, MassBal[t = 2:T*N_sim], S[t] == S[t-1] + s2hr*(q[t] - u[t])) # .- e*SA)
# (2) min max flow rate
@constraint(model, Release[t=1:T*N_sim], min_ut <= u[t] <= max_ut)
# (3) max solar capacity 
@constraint(model, SolarMax[t=1:T*N_sim], ps[t] == alpha_norm_w[t]*PS )
# (4) max feeder capacity 
@constraint(model, FeederMax[t=1:T*N_sim], ps[t] + ph[t] <= PF )
# (5) ramp rate constraint 
@constraint(model, RampRate[t=2:T*N_sim], RR_dn <= u[t] - u[t-1] <= RR_up)
# (6) daily water release contract
@constraint(model, Water_Contract, s2hr*sum(u) == Uw)
# (7) release to energy converstion 
@constraint(model, ReleaseEnergy[t=1:T*N_sim], ph[t] <= (eta*rho_w*g*s2hr*u[t]*(S[t]/SA))/ws2MWh)

# maximize revenue 
@objective(model, Max, R)

# initialize
# flattened for the whole week 
ps_s = zeros(T*N_sim, W_sim)
ph_s = zeros(T*N_sim, W_sim)
u_s = zeros(T*N_sim, W_sim)
S_s = zeros(T*N_sim, W_sim)
R_s = zeros(1, W_sim)
R_d = zeros(2, N_sim) # daily profit for both generators  

# Set the NonConvex parameter to 2
set_optimizer_attribute(model, "NonConvex", 2)

@printf("Optimization starts...\n")
@time begin 
for n = 1:W_sim # treat each week seperately 

optimize!(model)

global R_s[n] = value(R)      # objective_value(model);
global S_s[:,n] = value.(S)   # hourly
global ps_s[:,n] = value.(ps) # hourly
global ph_s[:,n] = value.(ph) # hourly
global u_s[:,n] = value.(u)   # hourly

# termination_status(model)
@printf("Finished Week %d, Cum Profit \$ %d Million, OptStatus: %s \n", n, sum(R_s)/1e6, termination_status(model))
@printf("Total FPV Profit \$ %d Million \n", sum(L.*ps_s)/1e6)
@printf("Total Hydropower Profit \$ %d Million \n \n", sum(L.*ph_s)/1e6)

@printf("Weekly Water Contract \$ %d m^3 \n \n", Uw)
@printf("Simulated Total Water Release \$ %d m^3 \n \n", sum(u_s)*s2hr)


for i = 1:N_sim
    s = (i - 1)*T + 1;
    f = i*T;

    Rfpv = sum(L[s:f].*ps_s[s:f]);   # revenue ($)
    Rhp = sum(L[s:f].*ph_s[s:f]);    # revenue ($)
    R_d[1,i] = Rfpv;                 # store values (1: FPV)
    R_d[2,i] = Rhp;                  # store values (2: HP)

    @printf("Day %d FPV Profit \$ %d \n", i, Rfpv)
    @printf("Day %d Hydropower Profit \$ %d \n", i, Rhp)
    @printf("\n")
end

end
end

# ---------- PLOTS -------- # 

# Create directory for this run 
stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS")
dir = "./plots/"
path = dir * stamp
mkdir(path)

# Generation Plots
hpcap = (eta*rho_w*g*s2hr*(u_s.*S_s)/SA)/ws2MWh # hydraulic head capacity on hydropower generation, convert to MWh
gen_plots(path, T*N_sim, ps_s, PS, alpha_norm_w, ph_s, hpcap, PF)

# Water release
release_plots(path, T*N_sim, u_s, min_ut, max_ut)

# Water release overlaid with LMP
release_plots_LMP(path, T*N_sim, u_s, min_ut, max_ut, L)

# ----- STATIC PLOTS ---- #
# Netflow 
num_days = Dates.daysinmonth(parse(Int, year), month_num)
flow_plot(path, num_days, netflow_2d[1:num_days])

# LMP
LMP_plot(path, T*N_sim, L)

# ----- NEXT STEPS ---- #
# TO DO: read in evap data 
# TO DO: initial conditions for ramp rate
# TO DO: after Wsim > 1, check L, Q, S0, and Uw are being updated to next day via update coef/RHS (initial conditions = final from prev week)
