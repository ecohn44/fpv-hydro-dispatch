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

#result = @benchmark begin
# elapsed_time = @elapsed begin

global s2hr = 3600  # seconds in an hour (delta t)
global min_ut = cfs_to_m3s(5000)    # min daily release limit [m3/s]
global max_ut = cfs_to_m3s(25000)   # max daily release limit [m3/s] 
global RR_dn = cfs_to_m3s(-2500) # down ramp rate limit [m3/s]
global RR_up = cfs_to_m3s(4000)  # up ramp rate limit [m3/s]
global PF = 1300   # max feeder capacity [MW] (3.3 GW) 
global PS = 1000   # max solar field capacity [MW] (1 GW) 
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
revenue = 0; 
thetas = 0;

Uw = sum(daily_s.release[1:N]) # monthly water contract
price = RTP_s.MW[1:T*N]
q = dailyflow_to_hourly(daily_s.inflow[1:N], T) # inflow
alpha_s = repeat(alpha[:,m], N) # solar radiation 


## ----------- RUN RELAXED ----------- ##
ps_t = []
ph_t = []
u_t = []

# Search bounds for DV
L = 0.0    #minimum(price)
R = 500.0  #maximum(price)
theta = (R + L)/2   
error = .01
i = 1
max_iter = 20 

while abs(R - L) > error

    theta = (R + L)/2 

    u, p_s, p_h, U_sim = run_sim_partialL(T, N, price, q, alpha_s, V0, theta);

    if U_sim > Uw # need to increase penalty to release less water, raise lower bounds
        L = theta
    end 
    if U_sim < Uw # decrease penalty to release more water, lower upper bounds
        R = theta
    end  

    # iterate 
    i = i + 1

    thetas = theta
    revenue = sum(price.*(p_s + p_h))
    ps_t = p_s
    ph_t = p_h 
    u_t = u
end

## SIMULATION RESULTS FOR TUNED THETA
# Total Revenue 
println(revenue)

# Price of Water
println(thetas)

# Total Water Release
println(sum(u_t)*s2hr)

# Water Contract
println(Uw)

#end

# Computation Time
#println("Median time in ms: $(median(result.times) / 1_000_000) ms")
#println("Elapsed time: $elapsed_time s")

# Download PS, PH, UT
combined_data = hcat(ps_t, ph_t, u_t)
df = DataFrame(combined_data, :auto)
CSV.write("output/weekly_behavior.csv", df)

