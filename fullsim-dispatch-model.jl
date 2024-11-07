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
using Statistics
using XLSX
include("plots.jl")
include("functions.jl")
include("dataload.jl")

global s2hr = 3600  # seconds in an hour (delta t)
global min_ut = cfs_to_m3s(5000)    # min daily release limit [m3/s]
global max_ut = cfs_to_m3s(25000)   # max daily release limit [m3/s] 
global RR_dn = cfs_to_m3s(-2500) # down ramp rate limit [m3/s]
global RR_up = cfs_to_m3s(4000)  # up ramp rate limit [m3/s]
# global PF        # max feeder capacity [MW] (nameplate 3.4 GW) 
global PS = 1000   # max solar field capacity [MW]  
global eta = .775     # efficiency of release-energy conversion
global rho_w = 1000 # density of water [kg/m^3]
global g = 9.8      # acceleration due to gravity [m/s^2]
global a = 14.9837;      # hydraulic head parameter 1 
global b = 0.1321;    # hydraulic head parameter 2 

monthly_overlayplots = false;
weeklyplots = false;
make_path = false;
baseline_sim = false; 
print = false; 
DV_plot = false;

# -----------------  DATA LOAD  ----------------- #
println("--- SIMULATION BEGIN ---")

feeder = [500, 1000, 1300, 2000, 3000, 4000]
years = ["22", "23"]
months = range(1, 12)
num_m = length(months)
num_y = length(years)
num_f = length(feeder)
T = 24 # hours (time steps)


# Load in 2022 - 2023 data & subset for month and year
daily, alpha, RTP = fullsim_dataload();

# Simulation Outputs
DVs = zeros(Float64, num_m, num_y, num_f) 
revenue = zeros(Float64, num_m, num_y, num_f)
release = zeros(Float64, num_m, num_y, num_f)
output_headers = []
input_headers = ["2022 LMP", "2023 LMP", "2022 WC", "2023 WC", "2022 Inflow", "2023 Inflow"]

# Simulation Inputs
avg_LMP = zeros(Float64, num_m, num_y)
water_contract = zeros(Float64, num_m, num_y) 
total_inflow = zeros(Float64, num_m, num_y)

for f in feeder
    global PF = f
    pf_ind = findfirst(x -> x == f, feeder)

    for y in years
        y_ind = findfirst(==(y), years)
        label = "20"*y*" PF"*string(f)
        push!(output_headers, label)

        for m in months

            # Monthly Datasets
            daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
            RTP_s = filter(row -> row[:Year] == y && row[:Month] == string(m), RTP)
            N = nrow(daily_s) # number of days in the month 
            alpha_s = repeat(alpha[:,m], N) # available solar radiation 
            V0 = daily_s.storage[1] # initial storage conditions
            q = dailyflow_to_hourly(daily_s.inflow, T) # inflow
            price = RTP_s.MW # price

            ## Outflow Water Contract
            Uw = sum(daily_s.release)

            println("Processing Month ", m)
            @printf("Water Contract %d \n", Uw)


            ## OPTIMIZATION
            # Run baseline multi-period simulation 
            if baseline_sim
                global min_Vt = V0 - T*N*s2hr*max_ut # min reservoir levels 
                u_b, ps_b, ph_b, f0_b, U_sim_b = run_sim(T, N, price, q, alpha_s, Uw, V0)
                @printf("Baseline Policy Revenue: \$ %d \n", f0_b)
                @printf("Water Release (Baseline): %d m3 \n", U_sim_b)
            end

            # Run partially relaxed formulation
            u, p_s, p_h, theta = bst_sim(T, N, price, q, alpha_s, V0, Uw)

            revenue[m,y_ind,pf_ind] = sum(price.*(p_s + p_h))
            release[m,y_ind,pf_ind] = sum(u)*s2hr
            DVs[m,y_ind,pf_ind] = theta

            avg_LMP[m,y_ind] = mean(price)
            total_inflow[m,y_ind] = sum(q*s2hr)
            water_contract[m,y_ind] = Uw
        end
    end
end

# ---------- SAVE CSV -------- #
## OUTPUT DATA
# save theta for 2022, 2023, PF
flat_DVs = reshape(DVs, size(DVs, 1), size(DVs, 2) * size(DVs, 3))
DV_df = DataFrame(flat_DVs, Symbol.(output_headers))
CSV.write("output/DVs.csv", DV_df)

# save revenue for 2022, 2023, PF
flat_rev = reshape(revenue, size(revenue, 1), size(revenue, 2) * size(revenue, 3))
rev_df = DataFrame(flat_rev, Symbol.(output_headers))
CSV.write("output/rev.csv", rev_df)

## INPUT DATA
# save LMP, inflow, outflow for 2022, 2023
combined_data = hcat(avg_LMP, water_contract, total_inflow)
input_df = DataFrame(combined_data, Symbol.(input_headers))
CSV.write("output/input.csv", input_df)


# ---------- SUMMARY -------- #
if print
    for m in months
        @printf("\nSummary for Month: %d \n", m)
        @printf("Relaxed Policy Revenue: \$ %d \n", revenue[m])
        @printf("Water Contract: %d m3 \n", water_contract[m])
        @printf("Water Release (Relaxed): %d m3 \n", release[m])
        @printf("Dual Value: %d \n", DVs[m])
    end
end

# ---------- PLOTS -------- #

if DV_plot
    scatter(DVs, label = "Theta (Dual Value)", legend = :outertopright)
    plot!(DVs, color=:blue, label = "")
    plot!(avg_LMP, label = "Avg LMP")
    title!("20"*y*" Water Contract Dual Values per Month")
    xaxis!("Month")
    yaxis!("Theta")
end

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

