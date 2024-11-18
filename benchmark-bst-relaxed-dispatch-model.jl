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
include("/Users/elizacohn/Desktop/fpv-hydro-dispatch/functions.jl")

# result = @benchmark begin

global s2hr = 3600  # seconds in an hour (delta t)
global min_ut = s2hr*cfs_to_m3s(5000)    # min daily release limit [m3/s] to [m3/hr] 
global max_ut = s2hr*cfs_to_m3s(25000)   # max daily release limit [m3/s] to [m3/hr] 
global RR_dn = s2hr*cfs_to_m3s(2500) # down ramp rate limit [m3/s] to [m3/hr] 
global RR_up = s2hr*cfs_to_m3s(4000)  # up ramp rate limit [m3/s] to [m3/hr] 
global P = 1300   # max feeder capacity [MW] (3.3 GW) 
global S = 1000   # max solar field capacity [MW] (1 GW) 
global eta = .775     # efficiency of release-energy conversion
global rho_w = 1000 # density of water [kg/m^3]
global g = 9.8      # acceleration due to gravity [m/s^2]
global a = 14.9837;      # hydraulic head parameter 1 
global b = 0.1321;     # hydraulic head parameter 2 

# -----------------  DATA LOAD  ----------------- #

y = "22"
m = 1
T = 24 
N = 7

# Load in 2022 - 2023 data
daily, alpha, RTP = fullsim_dataload();
daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
RTP_s = filter(row -> row[:Year] == y && row[:Month] == string(m), RTP)
V0 = daily_s.storage[1] # initial storage conditions
Uw = sum(daily_s.release[1:N]) # monthly water contract
price = RTP_s.MW[1:T*N]
q = dailyflow_to_hourly(daily_s.inflow[1:N], T) # inflow
alpha_s = repeat(alpha[:,m], N) # solar radiation 

## ----------- RUN RELAXED ----------- ##

# Search bounds for DV
L = 0.0    #minimum(price)
R = 10.0  #maximum(price)
theta = (R + L)/2   
error = .00000001
i = 1
max_iter = 50 

s = zeros(Float64, T*N) 
h = zeros(Float64, T*N) 
u = zeros(Float64, T*N) 
V = zeros(Float64, T*N) 
U_sim = 0
#theta = 0.03862637095153332


while abs(R - L) > error

    println(theta)
    theta = (R + L)/2 

    for t in range(1,T*N)
        
        # solar capacity 
        s[t] = min(alpha_s[t]*S, P);
        
        # set initial conditions
        if t == 1
            V_prev = V0
            u_prev = 0
        else
            u_prev = u[t-1]
            V_prev = V[t-1]
        end

        # phi = a*((V_prev)^b)
        # phi = a*((V0)^b)
        phi = 341.776
        theta_hat = (1/3.6e9)*price[t]*eta*g*rho_w*phi # ~0.039

        # water release
        if theta_hat > theta
            u_hat = (P-s[t])/((1/3.6e9)*eta*g*rho_w*phi) # ~1e7
        else
            u_hat = 0
        end

        u[t] = max(u_prev - RR_dn, min_ut, min(max_ut, u_prev + RR_up, u_hat))
        # mass balance and hydro generation
        V[t] = V_prev + q[t] - u[t]
        h[t] = min(P-s[t], eta*g*rho_w*u[t]*phi/3.6e9)
        #h[t] = min(P-s[t], eta*g*rho_w*u[t]*a*((V[t])^b)/3.6e9)
    end

    U_sim = sum(u)

    
    if U_sim > Uw # need to increase penalty to release less water, raise lower bounds
        L = theta
    end 
    if U_sim < Uw # decrease penalty to release more water, lower upper bounds
        R = theta
    end  

    # iterate 
    i = i + 1
# end


## SIMULATION RESULTS FOR TUNED THETA
# println("Simulation Results: ")
# Total Revenue 
revenue = sum(price.*(s + h))
println(revenue)

# Price of Water
println(theta)

# Total Water Release
println(U_sim)


end

# Computation Time
# println("Median time in ms: $(median(result.times) / 1_000_000) ms")

# Download PS, PH, UT
#combined_data = hcat(ps_t, ph_t, u_t)
#df = DataFrame(combined_data, :auto)
#CSV.write("output/relaxed_weekly_behavior.csv", df)

