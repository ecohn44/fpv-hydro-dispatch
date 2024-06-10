using JuMP
using Gurobi
using MAT
using Printf
using CSV
using DataFrames
using Dates
using Plots
using Base.Filesystem
using Ipopt


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
storage_path = string("data/USBR-lakemeadpowell-combinedvolume-",month,year,".csv");
hpcap_path = "data/hpcap.csv";

# Read in Local Marginal Price for Lake Mead
RTP = DataFrame(CSV.File(LMP_path));
RTP_2d = [row.LMP for row in eachrow(RTP)];

# Read in system inflow (lake powell inflow)
netflow = DataFrame(CSV.File(netflow_path));
# convert to m3/s from cfs
netflow_2d = [cfs_to_m3s(row["powell_inflow (cfs)"]) for row in eachrow(netflow)];

# Read in solar radiation capacity
alpha_raw = DataFrame(CSV.File(solarrad_path));
alpha_2d = [row[month] for row in eachrow(alpha_raw)];
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
storage_2d = [af_to_m3(row["total storage (ac-ft)"]) for row in eachrow(storage)];

# Read hpcap
hpcap = DataFrame(CSV.File(hpcap_path)).Column1;


# ----------------- REDUCTION FACTOR ----------------- #
# beta = .5

# ----------------- SIMULATION PARAMETERS ----------------- #
## Overview: Power contractors operate on weekly subcontracts. 
# Therefore, will let the simulation run every hour of the week before resetting for the next water contract 

T = 24; # time step per day
N_sim = 1; # days per week
W_sim = 1; # number of weeks 
s2hr = 3600  # seconds in an hour (delta T)

# Extract data for time period of interest
L = RTP_2d[1:T*N_sim,1];    # hourly LMP [$/MWh]
Q = netflow_2d[1:N_sim,1];  # daily system net flow into [m3]
U = release_2d[1:N_sim,1]   # daily total release contract [m3]
S = storage_2d[1:N_sim,1]   # historical daily storage levels [m3]
q = dailyflow_to_hourly(Q, T); # hourly system net flow (constant over each day) [m3]
alpha_norm_w = repeat(alpha_norm, N_sim) # replicate available solar capacity for every day of the sim. 
alpha = alpha_norm_w[1:T*N_sim]

V0 = S[1]   # initial reservoir conditions [~1.9 e10 m3]
e = 0    # evaporation constant [m/day]
min_ut = cfs_to_m3s(5000)    # min daily release limit [m3/s]
max_ut = cfs_to_m3s(25000)   # max daily release limit [m3/s] 
min_Vt = S0 - T*N_sim*s2hr*max_ut # min reservoir levels
RR_dn = cfs_to_m3s(-2500) # down ramp rate limit [m3/s]
RR_up = cfs_to_m3s(4000)  # up ramp rate limit [m3/s]
PF = 3300   # max feeder capacity [MW] (3.3 GW) 
PS = 1000   # max solar field capacity [MW] (1 GW) 
SA = 1.3e9   # surface area of both reservoirs [m^2] (504 sq miles)
eta = .9     # efficiency of release-energy conversion
rho_w = 1000 # density of water [kg/m^3]
g = 9.8      # acceleration due to gravity [m/s^2]

Uw = sum(U[1:N_sim]); # Weekly Water contract

# ------------ HYDRAULIC HEAD VARS ------------ #
a = 15;
b = 0.13; 

# ----------------- OPTIMIZATION  ----------------- #

# Create the optimization model
model = Model(Ipopt.Optimizer)

# Define variables
@variable(model, V[1:T])
set_lower_bound.(V, min_Vt)
@variable(model, p_h[1:T] >= 0)
@variable(model, p_s[1:T] >= 0)
@variable(model, u[1:T])
set_upper_bound.(u, max_ut)
set_lower_bound.(u, min_ut)

# Initial conditions
@constraint(model, V[1] == V0)

# Objective function
@objective(model, Max, sum(L .* (p_h + p_s)))

# Constraints
@constraint(model, [t in 2:T], V[t] == V[t-1] + s2hr*(q[t] - u[t]))
@constraint(model, s2hr*sum(u) == Uw)
@constraint(model, [t in 1:T], p_h[t] <= (eta * g * rho_w * u[t] * a * (V[t]^b))/1e6)
@constraint(model, [t in 1:T], min_ut <= u[t] <= max_ut)
@constraint(model, [t in 2:T], RR_dn <= u[t] - u[t-1] <= RR_up)
@constraint(model, [t in 1:T], 0 <= p_s[t] <= alpha[t]*PS)
@constraint(model, [t in 1:T], 0 <= p_s[t] + p_h[t] <= PF)

# Solve the optimization problem
optimize!(model)

# Print the results
println("Optimal V: ", value.(V))
println("Optimal p_h: ", value.(p_h))
println("Optimal p_s: ", value.(p_s))
println("Optimal u: ", value.(u))
println("Optimal objective value: ", objective_value(model))