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
include("functions.jl")



# ----------------- DATA LOAD / SIMULATION PARAMETERS  ----------------- #

# 1 WEEK SIMULATION FOR JAN 1-7 2022
year = "2022";
month = "January";

T = 24; # time step per day
N = 1; # days per week
s2hr = 3600  # seconds in an hour (delta t)

L, U, S, q, alpha = load_data(year, month, T, N)

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
a = 15;
b = 0.13; 


# ------------ MODEL PREDICTIVE CONTROL ------------ #

# Set moving horizon
K = 4

# Loop through total horizon, k is start index for reference 
for k = 1:(T*N - K)
    
    @printf("\n ///////////// Iteration: %d ///////////// \n" , k)

    # Determine day and index daily variables
    d = Int(floor(k/24) + 1)
    U_h = U[d]
    V0 = S[d]   # initial reservoir conditions [~1.9 e10 m3]

    # Index input vectors
    L_h = L[k:k+K-1] # perfect prediction
    q_h = q[k:k+K-1] # perfect prediction
    alpha_h = alpha[k:k+K-1] #perfect prediction

    Uw = sum(U_h)/6; # Weekly Water contract --> 4hr contract
    VT = V0 + sum(q_h) - Uw # Terminal volume conditions

    # T = 1, N = K for limited subhorizon
    run_sim(1, K, L_h, q_h, alpha_h, min_Vt, max_ut, min_ut, RR_up, RR_dn, PF, PS, V0, s2hr, eta, g, rho_w, a, b)
end


# ---------- PLOTS -------- # 

# Create directory for this run 
#stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS") ;
#dir = "./plots/" ;
#path = dir * stamp * " IPOPT";
#mkdir(path)

# Generation Plots
#hpcap = (eta * rho_w *g * a * value.(u) .* (value.(V).^b))/1e6 
#gen_plots(path, T*N, value.(p_s), PS, alpha_norm_w, value.(p_h), hpcap, PF)

# Water release
#release_plots(path, T*N, value.(u), min_ut, max_ut)

# Water release/generation overlaid with electicity price
#release_plots_LMP(path, T*N, value.(u), min_ut, max_ut, L)
#gen_plots_LMP(path, T*N, value.(p_s) + value.(p_h), PF, L)

