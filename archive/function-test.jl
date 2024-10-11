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

## TEMP PRICE LOAD 
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
daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
N = 31 #nrow(daily_s) 
price = RTP_2d[1:T*N,1];   
y_num = parse(Int, y) + 2000
alpha_s = repeat(alpha[:,m], N)
V0 = daily_s.storage[1] 
Uw = sum(daily_s.release[1:N]) 
q = dailyflow_to_hourly(daily_s.inflow, T)[1:T*N] 

#### --------- OPTIMIZATION ---------- ###
theta = 164 

u, p_s, p_h, V, f0, U_sim = run_sim_partialL(T, N, price, q, alpha_s, V0, theta);

## Plot ramp rate
plot(u[2:end]-u[1:end-1], label = "Simulated Ramp Rate", legend = :outertopright)   
hline!([RR_up], color=:red, linestyle=:dash, label="Max Ramp Rate")
hline!([RR_dn], color=:red, linestyle=:dash, label="Min Ramp Rate")
title!("Relaxed Ramp Rate for Simulated Policy")

