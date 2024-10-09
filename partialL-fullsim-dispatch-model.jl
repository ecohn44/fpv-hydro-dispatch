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


plots = false;
make_path = false;
search = true;

# ----------------- UNIT CONVERSION ----------------- #

# 1 WEEK SIMULATION FOR JAN 1-7 2022
year = "2022";
month = "January";

T = 24; # hours (constant, time steps)
N = 7; # days (variable)
s2hr = 3600  # seconds in an hour (delta t)

price, U, S, q, alpha = load_data(year, month, T, N)

V0 = S[1]   # initial reservoir conditions [~1.9 e10 m3]
Uw = sum(U[1:N]); # Water contract for N days
e = 0    # evaporation constant [m/day]
min_ut = cfs_to_m3s(5000)    # min daily release limit [m3/s]
max_ut = cfs_to_m3s(25000)   # max daily release limit [m3/s] 
RR_dn = cfs_to_m3s(-2500) # down ramp rate limit [m3/s]
RR_up = cfs_to_m3s(4000)  # up ramp rate limit [m3/s]
PF = 1200   # max feeder capacity [MW] (3.3 GW) 
PS = 1000   # max solar field capacity [MW] (1 GW) 
SA = 1.3e9   # surface area of both reservoirs [m^2] (504 sq miles)
eta = .9     # efficiency of release-energy conversion
rho_w = 1000 # density of water [kg/m^3]
g = 9.8      # acceleration due to gravity [m/s^2]
a = 15;      # hydraulic head parameter 1 
b = 0.13;    # hydraulic head parameter 2 

# ----------------- OPTIMIZATION  ----------------- #

# Search bounds for DV
eps = 2;
L = 0              # minimum(L) - eps;
R = 1000           # maximum(L) + eps;
theta = (R + L)/2   
error = 1
i = 1
max_iter = 10 # TO DO: make plots that subset on actual value of i

thetas = zeros(Float64, max_iter)
f0s = zeros(Float64, max_iter)
U_sims = zeros(Float64, max_iter)

if search 
    while R - L > error
        @printf("Iteration: %d \n", i)
        @printf("Upper Bound: %d \n", R)
        @printf("Lower Bound: %d \n", L)
        @printf("Theta: %d \n", theta)

        theta = (R + L)/2 
        u, p_s, p_h, V, f0, U_sim = run_sim_partialL(T, N, price, q, alpha, max_ut, min_ut, RR_up, RR_dn, PF, PS, V0, s2hr, eta, g, rho_w, a, b, theta)

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
end


# ---------- PLOTS -------- #

if make_path
    # Create directory for this run 
    stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS") ;
    dir = "./plots/" ;
    path = dir * stamp * " LPartial";
    mkdir(path)
end 

if plots
    # Convergence Plot
    iters_plot(path, max_iter, thetas, L"\theta_k")
    
    # Generation Plots
    hpcap = (eta * rho_w *g * a * value.(u) .* (value.(V).^b))/1e6 
    gen_plots(path, T*N, value.(p_s), PS, alpha_norm_w, value.(p_h), hpcap, PF)

    # Water release
    release_plots(path, T*N, value.(u), min_ut, max_ut)

    # Water release/generation overlaid with electicity price
    release_plots_LMP(path, T*N, value.(u), min_ut, max_ut, L)
    gen_plots_LMP(path, T*N, value.(p_s) + value.(p_h), PF, L)

    # Plot historically simulated dual values
    duals_plot(path, T*N-1, theta, L"\theta_t", "Mass Balance")
    duals_plot(path, T*N-1, theta_diff, L"\theta_t - \theta_{t-1}", "Mass Balance Moving Difference")
    duals_plot(path, T*N, mu10, L"\mu_{t,10}", "Energy from Water Release")

    # Plot Lagrangian policy
    policy_plot(path, T*N, policy)
    overlay_policy_plot(path, T*N, policy, value.(u))
    overlay_policy_plot_solar(path, T*N, policy, value.(p_s))

    head = a * (value.(V).^b) # [m]
    head_deriv = a *b  * value.(V).^(b-1)
    hhead_plots(path, T*N, head, head_deriv)
end