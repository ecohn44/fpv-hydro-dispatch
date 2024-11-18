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
using LaTeXStrings
using BenchmarkTools
include("/Users/elizacohn/Desktop/fpv-hydro-dispatch/dataload.jl")

global s2hr = 3600  # seconds in an hour (delta t)
global min_ut = s2hr*cfs_to_m3s(5000)    # min daily release limit [m3/s] to [m3/hr] 
global max_ut = s2hr*cfs_to_m3s(25000)   # max daily release limit [m3/s] to [m3/hr] 
global RR_dn = s2hr*cfs_to_m3s(-2500) # down ramp rate limit [m3/s] to [m3/hr] 
global RR_up = s2hr*cfs_to_m3s(4000)  # up ramp rate limit [m3/s] to [m3/hr] 
global PF = 1300   # max feeder capacity [MW] (3.3 GW) 
global PS = 1000   # max solar field capacity [MW] (1 GW) 
global SA = 1.3e9   # surface area of both reservoirs [m^2] (504 sq miles)
global eta = .775     # efficiency of release-energy conversion
global rho_w = 1000 # density of water [kg/m^3]
global g = 9.8      # acceleration due to gravity [m/s^2]
global a = 14.9837;      # hydraulic head parameter 1 
global b = 0.1321;     # hydraulic head parameter 2 

# -----------------  TIMING ----------------- #
obj = 0;
U_sim = 0;
DV = 0;

# result = @benchmark begin

# -----------------  DATA LOAD  ----------------- #

# years = ["22", "23"]
# months = range(1,12)
y = "22"
m = 1
T = 24 
N = 7

# Load in 2022 - 2023 data
daily, alpha, RTP = fullsim_dataload();
daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
RTP_s = filter(row -> row[:Year] == y && row[:Month] == string(m), RTP)

V0 = daily_s.storage[1] # initial storage conditions for month 
min_Vt = V0 - T*N*max_ut # min reservoir levels 
Uw = sum(daily_s.release[1:N]) # monthly water contract
price = RTP_s.MW[1:T*N] # time varying
q = dailyflow_to_hourly(daily_s.inflow[1:N], T) # inflow
alpha_s = repeat(alpha[:,m], N) # solar radiation 

## ----------- RUN BASELINE ----------- ##

# Create the optimization model
model = Model(Gurobi.Optimizer)
# model = Model(Ipopt.Optimizer)
set_silent(model) 

# Define variables
@variable(model, V[1:T*N])
set_lower_bound.(V, min_Vt)
@variable(model, p_h[1:T*N] >= 0)
@variable(model, p_s[1:T*N] >= 0)
@variable(model, u[1:T*N])
set_upper_bound.(u, max_ut)
set_lower_bound.(u, min_ut)

# Initial conditions
@constraint(model, MassBalInit, V[1] == V0)
@constraint(model, RampRateInit, u[1] == min_ut)

# Objective function
@objective(model, Max, sum(price .* (p_h + p_s)))

# Constraints
@constraint(model, MassBal[t in 2:T*N], V[t] == V[t-1] + (q[t] - u[t]))
@constraint(model, ReleaseEnergy[t in 1:T*N], p_h[t] <= (eta * g * rho_w * u[t] * 341.776)/3.6e9)
# @constraint(model, ReleaseEnergy[t in 1:T*N], p_h[t] <= (eta * g * rho_w * u[t] * a * (V[t]^b))/(3.6e9))
@constraint(model, Release[t in 2:T*N], min_ut <= u[t] <= max_ut)
@constraint(model, RampRate[t in 2:T*N], RR_dn <= u[t] - u[t-1] <= RR_up)
@constraint(model, SolarCap[t in 1:T*N], 0 <= 1000*p_s[t] <= 1000*alpha_s[t]*PS)
@constraint(model, FeederCap[t in 1:T*N], 0 <= p_s[t] + p_h[t] <= PF)
@constraint(model, WaterContract, sum(u) == Uw)  

# Solve the optimization problem
optimize!(model)

# Revenue
obj = objective_value(model);
println(obj)

# Total Water Release
U_sim = sum(value.(u));
println(U_sim)

# Water Contract
#println(Uw)

# end 

# Computation Time
# println("Median time in ms: $(median(result.times) / 1_000_000) ms")

# Download PS, PH, UT
#combined_data = hcat(ps_t, ph_t, u_t)
#df = DataFrame(combined_data, :auto)
#CSV.write("output/baseline_weekly_behavior.csv", df)

