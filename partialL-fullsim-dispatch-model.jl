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


weeklyplots = false;
make_path = false;
search = false;

# Temp
year = "2022";
month = "January";
month_num = 1;

LMP_path = string("data/LMP-meads-2-N101-",month,year,".csv");
RTP = DataFrame(CSV.File(LMP_path));
RTP_2d = [row.LMP for row in eachrow(RTP)];
price = RTP_2d[1:T*N,1];

# -----------------  DATA LOAD  ----------------- #

years = ["22"] #, "23"]
months = range(1,1) #12)
T = 24 # hours (time steps)

# Load in 2022 - 2023 data
daily, alpha, _ = fullsim_dataload()

for y in years

    for m in months
        # subset daily data
        daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
        N = nrow(daily_s) # number of days in the month 

        # subset price data
        y_num = parse(Int, y) + 2000
        # RTP_s = filter(row -> row[:Year] == y_num && row[:Month] == m, RTP)

        # subset and duplicate radiation data
        alpha_s = repeat(alpha[:,m], N)

        ##  OPTIMIZATION SETUP  
        V0 = daily_s.storage[1] # initial storage conditions
        Uw = sum(daily_s.release) # monthly water contract
        q = dailyflow_to_hourly(daily_s.inflow, T) # inflow
        # price = RTP_s.LMP # price 

        ## OPTIMIZATION
        # Run baseline multi-period simulation 
        u_b, ps_b, ph_b, f0_b, U_sim_b = run_sim(T, N, price, q, alpha_s, Uw, V0)

        # Run partially relaxed formulation
        u, p_s, p_h, f0, U_sim, theta, i = bst_sim(T, N, price, q, alpha_s, V0, Uw)
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