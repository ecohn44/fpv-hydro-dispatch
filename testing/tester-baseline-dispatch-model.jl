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
include("/Users/elizacohn/Desktop/fpv-hydro-dispatch/dataload.jl")

global s2hr = 3600  # seconds in an hour (delta t)
global min_ut = cfs_to_m3s(5000)    # min daily release limit [m3/s]
global max_ut = cfs_to_m3s(25000)   # max daily release limit [m3/s] 
global RR_dn = cfs_to_m3s(-2500) # down ramp rate limit [m3/s]
global RR_up = cfs_to_m3s(4000)  # up ramp rate limit [m3/s]
global PF = 1200   # max feeder capacity [MW] (3.3 GW) 
global PS = 1000   # max solar field capacity [MW] (1 GW) 
global SA = 1.3e9   # surface area of both reservoirs [m^2] (504 sq miles)
global eta = .9     # efficiency of release-energy conversion
global rho_w = 1000 # density of water [kg/m^3]
global g = 9.8      # acceleration due to gravity [m/s^2]
global a = 15;      # hydraulic head parameter 1 
global b = 0.13;    # hydraulic head parameter 2 

# -----------------  DATA LOAD  ----------------- #

# years = ["22", "23"]
# months = range(1,12)
y = "22"
m = 1
T = 24 
num_days = 31
N = 1

# Load in 2022 - 2023 data
daily, alpha, RTP = fullsim_dataload();
daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
RTP_s = filter(row -> row[:Year] == "20"*y && row[:Month] == string(m), RTP)

V0 = daily_s.storage[1] # initial storage conditions for month 
revenue = zeros(Float64, num_days) # revenue for each horizon length 

for N in 1:num_days
    #N = 7

    min_Vt = V0 - T*N*s2hr*max_ut # min reservoir levels 
    Uw = sum(daily_s.release[1:N]) # monthly water contract
    L = repeat([RTP_s.MW[1]], T*N) # keep price constant for each time step
    q = dailyflow_to_hourly(daily_s.inflow[1:N], T) # inflow
    alpha_s = repeat(alpha[:,m], N) # solar radiation 

    ## ----------- RUN BASELINE ----------- ##

    # Create the optimization model
    model = Model(Gurobi.Optimizer)
    #model = Model(Ipopt.Optimizer)
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
    @objective(model, Max, sum(L .* (p_h + p_s)))

    # Constraints
    @constraint(model, MassBal[t in 2:T*N], V[t] == V[t-1] + s2hr*(q[t] - u[t]))
    @constraint(model, ReleaseEnergy[t in 1:T*N], p_h[t] <= (eta * g * rho_w * u[t] * a * (V0^b))/1e6)
    #@constraint(model, ReleaseEnergy[t in 1:T*N], p_h[t] <= (eta * g * rho_w * u[t] * a * (V[t]^b))/1e6)
    @constraint(model, Release[t in 2:T*N], min_ut <= u[t] <= max_ut)
    @constraint(model, RampRate[t in 2:T*N], RR_dn <= u[t] - u[t-1] <= RR_up)
    @constraint(model, SolarCap[t in 1:T*N], 0 <= 1000*p_s[t] <= 1000*alpha_s[t]*PS)
    @constraint(model, FeederCap[t in 1:T*N], 0 <= p_s[t] + p_h[t] <= PF)
    #@constraint(model, WaterContract, s2hr*sum(u) == Uw)  

    # Solve the optimization problem
    optimize!(model)

    obj = objective_value(model);
    #print(obj)  

    revenue[N] = obj
end

CSV.write("testing/baseline_revenue_gurobi.csv", DataFrame(Revenue = revenue))
