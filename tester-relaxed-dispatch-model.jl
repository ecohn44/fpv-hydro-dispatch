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
include("plots.jl")
include("functions.jl")
include("dataload.jl")

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

y = "22"
m = 1
T = 24 
num_days = 31
N = 1

# Load in 2022 - 2023 data
daily, alpha, RTP = fullsim_dataload();

daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
RTP_s = filter(row -> row[:Year] == "20"*y && row[:Month] == string(m), RTP)

V0 = daily_s.storage[1] # initial storage conditions
min_Vt = V0 - T*N*s2hr*max_ut # min reservoir levels 
revenue = zeros(Float64, num_days) # revenue for each horizon length 
theta = 131

for N in 1:num_days
    Uw = sum(daily_s.release[1:N]) # monthly water contract
    L = repeat([RTP_s.MW[1]], T*N) # keep price constant for each time step
    q = dailyflow_to_hourly(daily_s.inflow[1:N], T) # inflow
    alpha_s = repeat(alpha[:,m], N) # solar radiation 


    ## ----------- RUN RELAXED ----------- ##
    # Create the optimization model
    model = Model(Ipopt.Optimizer)
    set_silent(model) 

    # Initialize empty vectors to store mass balance and control vars
    V_s = zeros(Float64, T*N) # Reservoir volume
    u_s = zeros(Float64, T*N)
    ps_s = zeros(Float64, T*N)
    ph_s = zeros(Float64, T*N)

    # Create the optimization model
    model = Model(Gurobi.Optimizer)
    #model = Model(Ipopt.Optimizer)
    set_silent(model) # no outputs

    # Define variables
    @variable(model, p_h, lower_bound = 0) # Hydro generation
    @variable(model, p_s, lower_bound = 0) # Solar generation
    @variable(model, u) # Water release 
    @variable(model, R) # Revenue

    # Objective Function
    @constraint(model, Rev, R == L[1]*p_s + L[1]*p_h) # - theta*u)

    # Constraints
    @constraint(model, ReleaseEnergy, p_h - ((eta * g * rho_w * a * (V0^b))/1e6)*u <= 0 )
    #@constraint(model, ReleaseEnergy, p_h - ((eta * g * rho_w * a * (V_s[1]^b))/1e6)*u <= 0 )
    @constraint(model, Release, min_ut <= u <= max_ut)
    @constraint(model, RampRateDn, -u <= -(RR_dn + u_s[1]))
    @constraint(model, RampRateUp, u <= RR_up + u_s[1])
    @constraint(model, SolarCap, 1000*p_s <= 1000*alpha_s[1]*PS)
    @constraint(model, FeederCap, 0 <= p_s + p_h <= PF)

    # maximize revenue 
    @objective(model, Max, R)

    for t = 1:T*N
        # Update electricity price coefficients (note sign flip; s_n_c moves terms to one side)
        set_normalized_coefficient(Rev, p_s, -L[t])
        set_normalized_coefficient(Rev, p_h, -L[t])
        set_normalized_rhs(SolarCap, 1000*alpha_s[t]*PS)

        if t == 1 # fix initial condition for ramp rate
            set_normalized_rhs(RampRateUp, min_ut)
            set_normalized_rhs(RampRateDn, -min_ut)
            #set_normalized_coefficient(ReleaseEnergy, u, -((eta * g * rho_w * a * (V0^b))/1e6))
        else
            set_normalized_rhs(RampRateUp, RR_up + u_s[t-1])
            set_normalized_rhs(RampRateDn, -(RR_dn + u_s[t-1]))
            #set_normalized_coefficient(ReleaseEnergy, u, -((eta * g * rho_w * a * (V_s[t-1]^b))/1e6))
        end

        # Solve the optimization problem
        optimize!(model)

        # Update mass balance + decision variables
        if t == 1
            V_s[t] = V0 + s2hr*(q[t] - value.(u)) # m3
        else
            V_s[t] = V_s[t-1] + s2hr*(q[t] - value.(u)) # m3
        end
        u_s[t] = value.(u) # m3/s
        ps_s[t] = value.(p_s)
        ph_s[t] = value.(p_h)
    end

    obj = sum(L.*(ps_s + ph_s));
    revenue[N] = obj
end 

CSV.write("relaxed_revenue_gurobi.csv", DataFrame(Revenue = revenue))