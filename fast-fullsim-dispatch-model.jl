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
global min_ut = s2hr*cfs_to_m3s(5000)    # min daily release limit [m3/s] to [m3/hr] 
global max_ut = s2hr*cfs_to_m3s(25000)   # max daily release limit [m3/s] to [m3/hr] 
global RR_dn = s2hr*cfs_to_m3s(2500) # down ramp rate limit [m3/s] to [m3/hr] 
global RR_up = s2hr*cfs_to_m3s(4000)  # up ramp rate limit [m3/s] to [m3/hr] 
# global PF        # max feeder capacity [MW] (nameplate 3.4 GW) 
global S = 1000   # max solar field capacity [MW]  
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

# Storage arrays for solar + hydro 
ST = Float64[]
HT = Float64[]

feeder = [1300] # [500, 1000, 1300, 2000, 3000, 4000]
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

for P in feeder
    pf_ind = findfirst(x -> x == P, feeder)

    for y in years
        y_ind = findfirst(==(y), years)
        label = "20"*y*" PF"*string(P)
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

            # Run partially relaxed formulation
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

                    phi = a*((V_prev)^b)
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
                    h[t] = min(P-s[t], eta*g*rho_w*u[t]*a*((V[t])^b)/3.6e9)
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
            end

            append!(ST, s)
            append!(HT, h)

            revenue[m,y_ind,pf_ind] = sum(price.*(s + h))
            release[m,y_ind,pf_ind] = sum(u)
            DVs[m,y_ind,pf_ind] = theta

            avg_LMP[m,y_ind] = mean(price)
            total_inflow[m,y_ind] = sum(q)
            water_contract[m,y_ind] = Uw
        end
    end
end

# ---------- SAVE CSV -------- #

# OUTPUT GENERATION
gen_headers = ["solar", "hydro"]
gen_data = hcat(ST, HT)
gen_df = DataFrame(gen_data, Symbol.(gen_headers))
CSV.write("output/fastgen.csv", gen_df)

#=
## OUTPUT DATA
# save theta for 2022, 2023, PF
flat_DVs = reshape(DVs, size(DVs, 1), size(DVs, 2) * size(DVs, 3))
DV_df = DataFrame(flat_DVs, Symbol.(output_headers))
CSV.write("output/fastDVs.csv", DV_df)

# save revenue for 2022, 2023, PF
flat_rev = reshape(revenue, size(revenue, 1), size(revenue, 2) * size(revenue, 3))
rev_df = DataFrame(flat_rev, Symbol.(output_headers))
CSV.write("output/fastrev.csv", rev_df)

## INPUT DATA
# save LMP, inflow, outflow for 2022, 2023
combined_data = hcat(avg_LMP, water_contract, total_inflow)
input_df = DataFrame(combined_data, Symbol.(input_headers))
CSV.write("output/fastinput.csv", input_df)
=#

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

println("--- SIMULATION END ---")

