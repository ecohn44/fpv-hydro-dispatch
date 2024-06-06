using JuMP
using Gurobi
using MAT
using Printf
using CSV
using DataFrames
using Dates
using Plots
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

# ----------------- DATA LOAD ----------------- #
# 24 HOUR SIMULATION FOR JAN 1 2022
year = "2022"
month = "January"

LMP_path = string("data/LMP-meads-2-N101-",month,year,".csv");
netflow_path = string("data/netflow-daily-",year,"-USBR.csv");
solarrad_path = "data/lake-mead-monthly-solar-rad-capacity.csv";
release_path = string("data/mead-totalrelease-",month,year,".csv");
storage_path = string("data/USBR-lakemead-daily-",month,year,".csv");

# Read in Local Marginal Price for Lake Mead
RTP = DataFrame(CSV.File(LMP_path));
RTP_2d = [row.LMP for row in eachrow(RTP)];

# Read in system net flow (lake powell inflow - lake mead outflow)
netflow = DataFrame(CSV.File(netflow_path));
# convert to m3/s from cfs
netflow_2d = [cfs_to_m3s(row.net) for row in eachrow(netflow)];

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

T = 24; # time step per day
N_sim = 1; # number of days 

# Extract data for time period of interest
L = RTP_2d[1:T*N_sim,1];    # hourly LMP [$/MWh]
Q = netflow_2d[1:N_sim,1];  # daily system net flow into [m3]
U = release_2d[1:N_sim,1]   # daily total release contract [m3]
S = storage_2d[1:N_sim,1]   # storage levels [m3]
q = (Q/T).*ones(Float64, 24)   # hourly system net flow (constant over day) [m3]

e = 0    # evaporation constant [m/day]
min_ut = cfs_to_m3s(5000)    # min daily release limit [m3/s]
max_ut = cfs_to_m3s(25000)   # min daily release limit [m3/s]
RR_dn = cfs_to_m3s(-2500) # down ramp rate limit [m3/s]
RR_up = cfs_to_m3s(4000)  # up ramp rate limit [m3/s]
PF = 99999   # max feeder capacity (~inf) 
PS = 10000   # max solar field capacity [MW] (10 GW) 
SA = 1.3e9   # surface area of both reservoirs [m^2] (504 sq miles)
eta = .9     # efficiency of release-energy conversion
rho_w = 1000 # density of water [kg/m^3]
g = 9.8      # acceleration due to gravity [m/s^2]
s2hr = 3600  # seconds in an hour 

# ----------------- VARIABLES ----------------- #

# Daily parameters 
Ud = U[1]; # Water contract
S0 = S[1]   # initial reservoir conditions [m3]

# initialize optimization model
model = Model(Gurobi.Optimizer)
set_silent(model) # no outputs

@variable(model, ph[1:T], lower_bound = 0)      # total hydro discharge power [Mw]
@variable(model, ps[1:T], lower_bound = 0)      # total solar discharge power [MW]
@variable(model, Sd, lower_bound = 0)           # EOD reservoir storage level [m3]
@variable(model, u[1:T], lower_bound = 0)       # water release [m3/s]
@variable(model, R)                             # market revenue [$]

# ----------------- OPTIMIZATION  ----------------- #

# operations revenue
@constraint(model, Revenue, R == sum(L.*(ph + ps)))
# reservoir mass balance equation
@constraint(model, MassBal, Sd == S0 + sum(s2hr.*(q - u))) # .- e*SA)
# min max flow rate
@constraint(model, Release[t=1:T], min_ut <= u[t] <= max_ut)
# max solar capacity 
@constraint(model, SolarMax[t=1:T], ps[t] <= alpha_norm[t]*PS )
# max feeder capacity 
@constraint(model, FeederMax[t=1:T], ps[t] + ph[t] <= PF )
# ramp rate constraint 
@constraint(model, RampRate[t = 2:T], RR_dn <= u[t] - u[t-1] <= RR_up)
# daily water release contract
@constraint(model, Water_Contract, s2hr*sum(u) <= Ud)
# release to energy converstion 
@constraint(model, ReleaseEnergy[t=1:T], ph[t] <= eta*rho_w*g*s2hr*u[t]*S0/SA)

# maximize revenue 
@objective(model, Max, R)

# initialize
ps_s = zeros(T, N_sim)
ph_s = zeros(T, N_sim)
u_s = zeros(T, N_sim)
S_s = zeros(1, N_sim)
R_s = zeros(1, N_sim)


@printf("Optimization starts...\n")
@time begin 
for n = 1:N_sim

optimize!(model)

global R_s[n] = value(R)      # objective_value(model);
global S_s[n] = value(Sd)     # daily 
global ps_s[:,n] = value.(ps) # hourly
global ph_s[:,n] = value.(ph) # hourly
global u_s[:,n] = value.(u)   # hourly

# termination_status(model)
@printf("Finished Day %d, Cum Profit \$ %d Million, OptStatus: %s \n", n, sum(R_s)/1e6, termination_status(model))
@printf("FPV Profit \$ %d Million \n", sum(L.*ps_s)/1e6)
@printf("Hydropower Profit \$ %d Million \n", sum(L.*ph_s)/1e6)


end
end

# ---------- PLOTS -------- # 
# eot_plots(T, ps_s, ph_s, u_s)

# Observations: 
# hydropower generation decoupled from water release