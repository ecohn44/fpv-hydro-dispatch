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

plots = false;

# ----------------- DATA LOAD / SIMULATION PARAMETERS  ----------------- #

# 1 WEEK SIMULATION FOR JAN 1-7 2022
year = "2022";
month = "January";

T = 16; # hours (time steps)
N = 1; # days per week
s2hr = 3600  # seconds in an hour (delta t)

L, U, S, q, alpha = load_data(year, month, T, N)

V0 = S[1]   # initial reservoir conditions [~1.9 e10 m3]
Uw = sum(U[1:N]); # Water contract for N days
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

# Create storage for optimal control 
u_star = [0.0 for n=1:T*N-1];
p_s_star =  [0.0 for n=1:T*N-1];
p_h_star = [0.0 for n=1:T*N-1];

# Loop through total horizon, k is start index for reference 
for k = 1:T*N-1

    # update horizon for end of period edge cases
    if k > T*N - K
        Kh = T*N-k
    else
        Kh = K
    end

    # prev value assignment
    if k == 1
        u_prev = 0
    else
        u_prev = u_star[k-1]
    end

    @printf("\n ///////////// Iteration: %d ///////////// \n" , k)

    # Determine day and index daily variables
    d = Int(floor(k/24) + 1)

    # Index input vectors
    L_h = L[k:k+Kh-1] # perfect prediction for price
    q_h = q[k:k+Kh-1] # perfect prediction for inflow 
    alpha_h = alpha[k:k+Kh-1] # historic available solar radiation 

    # T = 1, N = K for limited subhorizon
    u, p_s, p_h, mu11 = Lrun_sim(1, Kh, L_h, q_h, alpha_h, Uw, min_Vt, max_ut, min_ut, RR_up, RR_dn, PF, PS, V0, s2hr, eta, g, rho_w, a, b, k, u_prev, u_star)

    # Store control decision for t = 1
    u_star[k] = u[1]
    p_h_star[k] = p_h[1]
    p_s_star[k] = p_s[1]

    # testing 
    #if k == 2
    #    break
    #end

end


# ---------- PLOTS -------- # 

if plots
    # Create directory for this run 
    stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS") ;
    dir = "./plots/" ;
    path = dir * stamp * " IPOPT";
    mkdir(path)

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