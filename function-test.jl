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
include("functions.jl")

year = "2022";
month = "January";
month_num = 1;

LMP_path = string("data/LMP-meads-2-N101-",month,year,".csv");
RTP = DataFrame(CSV.File(LMP_path));
RTP_2d = [row.LMP for row in eachrow(RTP)];

# ----------------- DATA LOAD / SIMULATION PARAMETERS ----------------- #
y = "22"
m = 1

daily, alpha, _ = fullsim_dataload()

# subset daily data
daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
N = 31 #nrow(daily_s) # number of days in the month 
price = RTP_2d[1:T*N,1];    # hourly LMP [$/MWh]

# subset price data
y_num = parse(Int, y) + 2000
#RTP_s = filter(row -> row[:Year] == y_num && row[:Month] == m, RTP)

# subset and duplicate radiation data
alpha_s = repeat(alpha[:,m], N)

##  OPTIMIZATION SETUP  
V0 = daily_s.storage[1] # initial storage conditions
Uw = sum(daily_s.release[1:N]) # monthly water contract
q = dailyflow_to_hourly(daily_s.inflow, T)[1:T*N] # inflow
# price = RTP_s.LMP[1:T*N] # price 

# Search bounds for DV
eps = 2;
L = 0   #minimum(price)
R = 500  #maximum(price)
theta = (R + L)/2   
error = 1
i = 1
max_iter = 20 

thetas = zeros(Float64, max_iter)
f0s = zeros(Float64, max_iter)
U_sims = zeros(Float64, max_iter)

while abs(R - L) > error
    @printf("Iteration: %d \n", i)
    @printf("Upper Bound: %d \n", R)
    @printf("Lower Bound: %d \n", L)
    @printf("Theta: %d \n", theta)

    theta = (R + L)/2 

    u, p_s, p_h, V, f0, U_sim = run_sim_partialL(T, N, price, q, alpha_s, V0, theta)

    if U_sim > Uw # need to increase penalty to release less water, raise lower bounds
        L = theta
    end 
    if U_sim < Uw # decrease penalty to release more water, lower upper bounds
        R = theta
    end  

    # store  value of theta 
    thetas[i] = theta
    # store value of objective function 
    U_sims[i] = U_sim

    # iterate 
    i = i + 1
end

