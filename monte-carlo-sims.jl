using JuMP
using Gurobi
using MAT
using Printf
using CSV
using DataFrames
using Dates
using Plots
using Base.Filesystem
using LaTeXStrings
using Statistics
using XLSX
using StatsPlots
using Random
include("plots.jl")
include("dataload.jl")
include("monte-carlo-functions.jl")

global s2hr = 3600  # seconds in an hour (delta t)
global min_ut = s2hr*cfs_to_m3s(5000)    # min daily release limit [m3/s] to [m3/hr] 
global max_ut = s2hr*cfs_to_m3s(25000)   # max daily release limit [m3/s] to [m3/hr] 
global RR_dn = s2hr*cfs_to_m3s(2500) # down ramp rate limit [m3/s] to [m3/hr] 
global RR_up = s2hr*cfs_to_m3s(4000)  # up ramp rate limit [m3/s] to [m3/hr] 
global P = 1300       # max feeder capacity [MW] (nameplate 3.4 GW) 
global S = 1000   # max solar field capacity [MW]  
global eta = .775     # efficiency of release-energy conversion
global rho_w = 1000 # density of water [kg/m^3]
global g = 9.8      # acceleration due to gravity [m/s^2]
global a = 14.9837;      # hydraulic head parameter 1 
global b = 0.1321;    # hydraulic head parameter 2 


## Data Load 
println("--- SIMULATION BEGIN ---")
# Load in 2022 - 2023 data & subset for month and year
daily, alpha, RTP = fullsim_dataload();
mapes = [0 , 0.05, 0.1, 0.15, 0.2] # noise level parameter
phi_mc = 0.8 # noise hyperparameter
num_iters = 100 # Monte carlo simulations 
years = ["22", "23"]
months = range(1, 12)
num_m = length(months)
num_y = length(years)
T = 24 

total_revenue = zeros(Float64, num_iters, length(mapes))
total_release = zeros(Float64, num_iters, length(mapes))
avg_revenue = zeros(Float64, num_iters, length(mapes))
avg_release = zeros(Float64, num_iters, length(mapes))

# Load in water price from fast-fullsim-dispatch-model 
DVs_path = string("output/fastDVs.csv");
DVs = DataFrame(CSV.File(DVs_path));

# Select water price for designated feeder capacity
selected_cols = filter(name -> endswith(string(name), string(P)), names(DVs))
waterprice = hcat(DVs[:, selected_cols])

# Load in pre-tuned AR(1) noise sigmas from monte-carlo-tuning
sigma_path = string("data/monte-carlo-sigmas.csv")
sigmas = DataFrame(CSV.File(sigma_path));

## Set noise parameters 
for mae in mapes
    # Set random seed 
    Random.seed!(42)
    println("Running for target MAE = ", mae)
    mae_index = findfirst(==(mae), mapes)

    # Extract pre-tuned sigmas for MAE level 
    col_symbol = Symbol(mae)
    sigma_mae = col_symbol in propertynames(sigmas) ? sigmas[:, col_symbol] : zeros(3)

    ## Run Monte Carlo simulations 
    for i in 1:num_iters

        # Add noise to price 
        noisy_price = add_noise_norm(sigma_mae[1], phi_mc, length(RTP.MW), RTP.MW)
        RTP.noisy_price = noisy_price

        # Add noise to inflow
        noisy_inflow = add_noise_norm(sigma_mae[2], phi_mc, length(daily.inflow), daily.inflow)
        daily.noisy_inflow = noisy_inflow

        # Simulation Outputs
        water_contract = zeros(Float64, num_m, num_y) 
        revenue = zeros(Float64, num_m, num_y)
        release = zeros(Float64, num_m, num_y)

        for y in years
            y_ind = findfirst(==(y), years)

            for m in months

                # Monthly Datasets
                daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)
                RTP_s = filter(row -> row[:Year] == y && row[:Month] == string(m), RTP)
                N = nrow(daily_s) # number of days in the month 
                alpha_s = repeat(add_noise(sigma_mae[3], phi_mc, T, alpha[:,m]), N) # available noisy solar radiation 
                V0 = daily_s.storage[1] # initial storage conditions
                q = dailyflow_to_hourly(daily_s.noisy_inflow, T) # noisy inflow
                price = RTP_s.noisy_price # noisy price

                ## Outflow Water Contracts
                Uw = sum(daily_s.release)
                ## Decision Variables
                s = zeros(Float64, T*N) 
                h = zeros(Float64, T*N) 
                u = zeros(Float64, T*N) 
                V = zeros(Float64, T*N) 

                # Use monthly optimal price of water 
                theta = waterprice[m,y_ind]

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

                    phi = a*((V_prev)^b)
                    theta_hat = (1/3.6e9)*price[t]*eta*g*rho_w*phi # ~0.039

                    # water release
                    if theta_hat > theta
                        u_hat = (P-s[t])/((1/3.6e9)*eta*g*rho_w*phi) # ~1e7
                    else
                        u_hat = 0
                    end

                    # Stochastic feasability (max)
                    if (T*N - t)*max_ut <= Uw - sum(u) 
                        # Set release to max 
                        u_hat = max_ut
                    end  

                    u[t] = max(u_prev - RR_dn, min_ut, min(max_ut, u_prev + RR_up, u_hat))
                    
                    # Stochastic feasability (min)
                    if Uw - sum(u) <= 0 
                        # No water remaining
                        u[t] = 0
                    end  

                    # mass balance and hydro generation
                    V[t] = V_prev + q[t] - u[t]
                    h[t] = min(P-s[t], eta*g*rho_w*u[t]*a*((V[t])^b)/3.6e9)
                end

                revenue[m,y_ind] = sum(price.*(s + h))
                release[m,y_ind] = sum(u)
                water_contract[m,y_ind] = Uw
            end
        end

        # Sum total revenue 
        total_revenue[i,mae_index] = sum(revenue)
        total_release[i,mae_index] = sum(release)
        avg_revenue[i,mae_index] = mean(revenue)
        avg_release[i,mae_index] = mean(release)
    end
end
println("--- SIMULATION END ---")

boxplot(1:length(mapes), eachcol(total_revenue/1e6), legend=false, xlabel="MAPE Level (%)", ylabel="Total Revenue (\$ M)", title="Total Revenue vs MAPE")
boxplot(1:length(mapes), eachcol(total_release/1e6), legend=false, xlabel="MAPE Level (%)", ylabel="Total Release (M m3)", title="Total Release vs MAPE")

boxplot(1:length(mapes), eachcol(avg_revenue/1e6), legend=false, xlabel="MAPE Level (%)", ylabel="Avg Revenue (\$ M)", title="Avg Revenue vs MAPE")
boxplot(1:length(mapes), eachcol(avg_release/1e6), legend=false, xlabel="MAPE Level (%)", ylabel="Avg Release (M m3)", title="Avg Release vs MAPE")

CSV.write("output/monte-carlo-revenue-boxplot.csv", DataFrame(total_revenue, :auto))
CSV.write("output/monte-carlo-release-boxplot.csv", DataFrame(total_release, :auto))