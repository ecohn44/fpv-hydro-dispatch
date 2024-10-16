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

monthly_overlayplots = false;
weeklyplots = false;
make_path = false;

# -----------------  DATA LOAD  ----------------- #
println("--- SIMULATION BEGIN ---")

# years = ["22", "23"]
# months = range(1,12)
y = "22"
m = 1
T = 24 # hours (time steps)

# Load in 2022 - 2023 data & subset for month and year
daily, alpha, RTP = fullsim_dataload();
daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
RTP_s = filter(row -> row[:Year] == "20"*y && row[:Month] == string(m), RTP)
N = nrow(daily_s) # number of days in the month 
alpha_s = repeat(alpha[:,m], N)

##  OPTIMIZATION SETUP  
V0 = daily_s.storage[1] # initial storage conditions
global min_Vt = V0 - T*N*s2hr*max_ut # min reservoir levels 
Uw = sum(daily_s.release) # monthly water contract
q = dailyflow_to_hourly(daily_s.inflow, T) # inflow
price = RTP_s.MW # price 

## OPTIMIZATION
# Run baseline multi-period simulation 
u_b, ps_b, ph_b, f0_b, U_sim_b = run_sim(T, N, price, q, alpha_s, Uw, V0)

# Run partially relaxed formulation
u, p_s, p_h, theta = bst_sim(T, N, price, q, alpha_s, V0, Uw)

## MONTHLY SUMMARY
@printf("Summary for Month: %d \n", m)
@printf("Baseline Policy Revenue: \$ %d \n", f0_b)
@printf("Relaxed Policy Revenue: \$ %d \n", sum(price.*(p_s + p_h)))
@printf("Water Contract: %d m3 \n", Uw)
@printf("Water Release (Baseline): %d m3 \n", U_sim_b)
@printf("Water Release (Relaxed): %d m3 \n", sum(u)*s2hr)
@printf("Dual Value: %d \n", theta)

# ---------- PLOTS -------- #

if make_path
    # Create directory for this run 
    stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS") ;
    dir = "./plots/" ;
    path = dir * stamp * " LPartial";
    mkdir(path)
end 

if monthly_overlayplots
    # plot price of electricity 
    LMP_plot(path, T*N, price) 
    
    # plot overlay of baseline water release with price
    release_plots_LMP(path, T*N, u_b, min_ut, max_ut, price) 

    # plot overlay of water releases 
    release_overlay_plots(path, T*N, u_b, u, min_ut, max_ut)

    # plot overlay of solar generation 
    generation_overlay_plots(path, T*N, ps_b, p_s, "Solar", PS)

    # plot overlay of hydropower generation
    generation_overlay_plots(path, T*N, ph_b, p_h, "Hydropower", PF)

    # plot overlay of total generation
    generation_overlay_plots(path, T*N, ph_b + ps_b, p_h + p_s, "Joint", PF)
end

if weeklyplots
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

println("--- SIMULATION END ---")