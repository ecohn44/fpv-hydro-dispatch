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
include("/Users/elizacohn/Desktop/fpv-hydro-dispatch/functions.jl")

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
revenue = zeros(Float64, num_days) # revenue for each horizon length 
thetas = zeros(Float64, num_days)

for N in 1:num_days

    Uw = sum(daily_s.release[1:N]) # monthly water contract
    price = RTP_s.MW[1:T*N]
    q = dailyflow_to_hourly(daily_s.inflow[1:N], T) # inflow
    alpha_s = repeat(alpha[:,m], N) # solar radiation 


    ## ----------- RUN RELAXED ----------- ##

    # Search bounds for DV
    L = 0.0    #minimum(price)
    R = 500.0  #maximum(price)
    theta = (R + L)/2   
    error = 1 #.01
    i = 1
    max_iter = 20 

    while abs(R - L) > error

        theta = (R + L)/2 

        u, p_s, p_h, V, U_sim = run_sim_partialL(T, N, price, q, alpha_s, V0, theta);

        if U_sim > Uw # need to increase penalty to release less water, raise lower bounds
            L = theta
        end 
        if U_sim < Uw # decrease penalty to release more water, lower upper bounds
            R = theta
        end  

        # iterate 
        i = i + 1

        thetas[N] = theta
        revenue[N] = sum(price.*(p_s + p_h))
    end
        
end 

CSV.write("testing/relaxed_revenue_gurobi_optwc.csv", DataFrame(Revenue = revenue))
