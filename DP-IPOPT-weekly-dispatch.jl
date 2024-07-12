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
include("plots.jl")
include("DP-functions.jl")

# ----------------- DATA LOAD / SIMULATION PARAMETERS ----------------- #

# 1 WEEK SIMULATION FOR JAN 1-7 2022
year = "2022";
month = "January";

T = 24; # time step per day
N = 7; # days per week
s2hr = 3600  # seconds in an hour (delta t)

L, U, S, q, alpha = load_data(year, month, T, N)

Uw = sum(U[1:N]); # Weekly Water contract
V0 = S[1]   # initial reservoir conditions [~1.9 e10 m3]
e = 0    # evaporation constant [m/day]
min_ut = cfs_to_m3s(5000)    # min daily release limit [m3/s]
max_ut = cfs_to_m3s(25000)   # max daily release limit [m3/s] 
min_Vt = V0 - T*N*s2hr*max_ut # min reservoir levels
RR_dn = cfs_to_m3s(-2500) # down ramp rate limit [m3/s]
RR_up = cfs_to_m3s(4000)  # up ramp rate limit [m3/s]
PF = 1200   # max feeder capacity [MW] (3.3 GW) 
PS = 1000   # max solar field capacity [MW] (1 GW) 
SA = 1.3e9   # surface area of both reservoirs [m^2] (504 sq miles)
eta = .9     # efficiency of release-energy conversion
rho_w = 1000 # density of water [kg/m^3]
g = 9.8      # acceleration due to gravity [m/s^2]
a = 15; # HYDRAULIC HEAD VARS
b = 0.13; #  HYDRAULIC HEAD VARS

# ------------ RUN OPTIMIZATION ------------ #

run_sim(T, N, L, q, alpha, min_Vt, max_ut, min_ut, RR_up, RR_dn, PF, PS, V0, s2hr, eta, g, rho_w, a, b)