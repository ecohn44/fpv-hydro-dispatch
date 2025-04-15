using JuMP
using Gurobi

# decomposed
function run_sim_partialL(T, N, L, q, alpha, V0, theta)

    print = false; 

    # Initialize empty vectors to store mass balance and control vars
    V_s = zeros(Float64, T*N) # Reservoir volume
    u_s = zeros(Float64, T*N)
    ps_s = zeros(Float64, T*N)
    ph_s = zeros(Float64, T*N)

    # Create the optimization model
    model = Model(Gurobi.Optimizer)
    #model = Model(Ipopt.Optimizer)
    set_silent(model) # no outputs

    # Define variables
    @variable(model, p_h, lower_bound = 0) # Hydro generation
    @variable(model, p_s, lower_bound = 0) # Solar generation
    @variable(model, u) # Water release 
    @variable(model, R) # Revenue
    
    # Objective Function
    @constraint(model, Rev, R == L[1]*p_s + L[1]*p_h - theta*u)

    # Constraints
    @constraint(model, ReleaseEnergy, p_h - ((eta * g * rho_w * a * (V_s[1]^b))/1e6)*u <= 0 )
    @constraint(model, Release, min_ut <= u <= max_ut)
    @constraint(model, RampRateDn, -u <= -(RR_dn + u_s[1]))
    @constraint(model, RampRateUp, u <= RR_up + u_s[1])
    @constraint(model, SolarCap, 1000*p_s <= 1000*alpha[1]*PS)
    @constraint(model, FeederCap, 0 <= p_s + p_h <= PF)

    # maximize revenue 
    @objective(model, Max, R)

    for t = 1:T*N
        # Update electricity price coefficients (note sign flip; s_n_c moves terms to one side)
        set_normalized_coefficient(Rev, p_s, -L[t])
        set_normalized_coefficient(Rev, p_h, -L[t])
        set_normalized_rhs(SolarCap, 1000*alpha[t]*PS)

        if t == 1 # fix initial condition for ramp rate
            set_normalized_rhs(RampRateUp, min_ut)
            set_normalized_rhs(RampRateDn, -min_ut)
            set_normalized_coefficient(ReleaseEnergy, u, -((eta * g * rho_w * a * (V0^b))/1e6))
        else
            set_normalized_rhs(RampRateUp, RR_up + u_s[t-1])
            set_normalized_rhs(RampRateDn, -(RR_dn + u_s[t-1]))
            set_normalized_coefficient(ReleaseEnergy, u, -((eta * g * rho_w * a * (V_s[t-1]^b))/1e6))
        end

        # Solve the optimization problem
        optimize!(model)

        # Update mass balance + decision variables
        if t == 1
            V_s[t] = V0 + s2hr*(q[t] - value.(u)) # m3
        else
            V_s[t] = V_s[t-1] + s2hr*(q[t] - value.(u)) # m3
        end
        u_s[t] = value.(u) # m3/s
        ps_s[t] = value.(p_s)
        ph_s[t] = value.(p_h)
    end

    total_release = sum(u_s*s2hr)
    
    # Print the final results 
    if print
        @printf("Simulated Release: %d m^3 \n \n", total_release)
    end 

    # return optimal control vars (x3), volume state var, total revenue, total release
    return u_s, ps_s, ph_s, total_release
    
end

# baseline
function run_sim(T, N, L, q, alpha, Uw, V0)

    print = false; 

    # Create the optimization model
    model = Model(Gurobi.Optimizer)
    set_silent(model) # no outputs

    # Define variables
    @variable(model, V[1:T*N])
    set_lower_bound.(V, min_Vt)
    @variable(model, p_h[1:T*N] >= 0)
    @variable(model, p_s[1:T*N] >= 0)
    @variable(model, u[1:T*N])
    set_upper_bound.(u, max_ut)
    set_lower_bound.(u, min_ut)

    # Initial conditions
    @constraint(model, MassBalInit, V[1] == V0)
    @constraint(model, RampRateInit, u[1] == min_ut)
    
    # Objective function
    @objective(model, Max, sum(L .* (p_h + p_s)))

    # Constraints
    @constraint(model, MassBal[t in 2:T*N], V[t] == V[t-1] + s2hr*(q[t] - u[t]))
    @constraint(model, ReleaseEnergy[t in 1:T*N], p_h[t] <= (eta * g * rho_w * u[t] * a * (V[t]^b))/3.6e9)
    @constraint(model, Release[t in 2:T*N], min_ut <= u[t] <= max_ut)
    @constraint(model, RampRate[t in 2:T*N], RR_dn <= u[t] - u[t-1] <= RR_up)
    @constraint(model, SolarCap[t in 1:T*N], 0 <= 1000*p_s[t] <= 1000*alpha[t]*PS)
    @constraint(model, FeederCap[t in 1:T*N], 0 <= p_s[t] + p_h[t] <= PF)
    @constraint(model, WaterContract, s2hr*sum(u) == Uw)  

    # Solve the optimization problem
    optimize!(model)

    obj = objective_value(model);
    total_release = sum(value.(u))*s2hr

    if print 
        # Print the results
        @printf("Total Profit: \$ %d \n", obj)
        @printf("Total FPV Profit: \$ %d \n", sum(L.*value.(p_s)))
        @printf("Total Hydropower Profit: \$ %d \n \n", sum(L.*value.(p_h)))

        @printf("Weekly Water Contract: %d m^3 \n \n", Uw)
        @printf("Simulated Total Water Release: %d m^3 \n \n", total_release)

        # ---------- DUAL VALUES -------- # 
        println("Dual Values")
        println("Water Contract: ", -dual.(WaterContract)*s2hr)
    end 

    # return optimal control vars, revenue, and total release
    return value.(u), value.(p_s), value.(p_h), obj, total_release

end

function bst_sim(T, N, price, q, alpha_s, V0, Uw)
    print = false; 

    # Search bounds for DV
    L = 0    #minimum(price)
    R = 1200 #maximum(price)
    theta = (R + L)/2   
    error = .01
    i = 1

    ut_sim = []
    ps_sim = []
    ph_sim = []

    while abs(R - L) > error
        if print
            @printf("Iteration: %d \n", i)
            @printf("Upper Bound: %d \n", R)
            @printf("Lower Bound: %d \n", L)
            @printf("Theta: %d \n", theta)
        end

        theta = (R + L)/2 

        u, p_s, p_h, U_sim = run_sim_partialL(T, N, price, q, alpha_s, V0, theta);

        if U_sim > Uw # need to increase penalty to release less water, raise lower bounds
            L = theta
        end 
        if U_sim < Uw # decrease penalty to release more water, lower upper bounds
            R = theta
        end  

        # store control decisions 
        ut_sim = u
        ps_sim = p_s
        ph_sim = p_h

        # iterate 
        i = i + 1
    end

    # return optimal trajectory, optimal DV, and iterations to convergency 
    return ut_sim, ps_sim, ph_sim, theta
end