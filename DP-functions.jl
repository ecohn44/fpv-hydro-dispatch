function cfs_to_m3s(cfs)
    # Conversion factor: 1 cfs is approximately 0.0283168 m³/s
    conversion_factor = 0.0283168
    return cfs * conversion_factor
end

function af_to_m3(acft)
    # Conversion factor: 1 acre-foot is approximately 1233.48 cubic meters
    conversion_factor = 1233.48
    return acft * conversion_factor
end

function dailyflow_to_hourly(q, T)
    # Input: Vector q of average daily flow values
    # Output: Hourly flow vector held constant over the day
    hourly_q = []
    for i in q
        rpt = fill(i, T)
        append!(hourly_q, rpt)
    end
    return hourly_q
end

function run_sim(T, N, L, q, alpha, min_Vt, max_ut, min_ut, RR_up, RR_dn, PF, PS, V0, s2hr, eta, g, rho_w, a, b)

    # Create the optimization model
    model = Model(Ipopt.Optimizer)
    # model = Model(Gurobi.Optimizer)

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

    # End conditions
    # @constraint(model, WaterContract, V[T*N] >= V0 - Uw + s2hr*sum(q))

    # Objective function
    @objective(model, Max, sum(L .* (p_h + p_s)))

    # Constraints
    @constraint(model, MassBal[t in 2:T*N], V[t] == V[t-1] + s2hr*(q[t] - u[t]))
    @constraint(model, WaterContract, s2hr*sum(u) == Uw)  
    @constraint(model, ReleaseEnergy[t in 1:T*N], p_h[t] <= (eta * g * rho_w * u[t] * a * (V[t]^b))/1e6)
    @constraint(model, Release[t in 1:T*N], min_ut <= u[t] <= max_ut)
    @constraint(model, RampRate[t in 2:T*N], RR_dn <= u[t] - u[t-1] <= RR_up)
    @constraint(model, SolarCap[t in 1:T*N], 0 <= p_s[t] <= alpha[t]*PS)
    @constraint(model, FeederCap[t in 1:T*N], 0 <= p_s[t] + p_h[t] <= PF)

    # Solve the optimization problem
    optimize!(model)

    obj = objective_value(model);

    # Print the results
    @printf("Total Profit: \$ %d \n", obj)
    @printf("Total FPV Profit: \$ %d \n", sum(L.*value.(p_s)))
    @printf("Total Hydropower Profit: \$ %d \n \n", sum(L.*value.(p_h)))

    @printf("Weekly Water Contract: %d m^3 \n \n", Uw)
    @printf("Simulated Total Water Release: %d m^3 \n \n", sum(value.(u))*s2hr)
end

function load_data(year, month, T, N)
    LMP_path = string("data/LMP-meads-2-N101-",month,year,".csv");
    netflow_path = string("data/netflow-daily-",year,"-USBR.csv");
    solarrad_path = "data/lake-mead-monthly-solar-rad-capacity.csv";
    release_path = string("data/mead-totalrelease-",month,year,".csv");
    storage_path = string("data/USBR-lakemeadpowell-combinedvolume-",month,year,".csv");

    # Read in Local Marginal Price for Lake Mead
    RTP = DataFrame(CSV.File(LMP_path));
    RTP_2d = [row.LMP for row in eachrow(RTP)];

    # Read in system inflow (lake powell inflow)
    netflow = DataFrame(CSV.File(netflow_path));
    # convert to m3/s from cfs
    netflow_2d = [cfs_to_m3s(row["powell_inflow (cfs)"]) for row in eachrow(netflow)];

    # Read in solar radiation capacity
    alpha_raw = DataFrame(CSV.File(solarrad_path));
    alpha_2d = [row[month] for row in eachrow(alpha_raw)];
    min_val = minimum(alpha_2d);
    max_val = maximum(alpha_2d);
    alpha_norm = (alpha_2d .- min_val) / (max_val - min_val); # Normalize data 

    # Read in Water Contract Data
    release = DataFrame(CSV.File(release_path));
    # convert from ac-ft to m3
    release_2d = [af_to_m3(row.Result) for row in eachrow(release)];

    # Read in historical storage levels 
    storage = DataFrame(CSV.File(storage_path));
    # convert from ac-ft to m3
    storage_2d = [af_to_m3(row["total storage (ac-ft)"]) for row in eachrow(storage)];

    # Extract data for time period of interest
    L = RTP_2d[1:T*N,1];    # hourly LMP [$/MWh]
    Q = netflow_2d[1:N,1];  # daily system net flow into [m3]
    U = release_2d[1:N,1]   # daily total release contract [m3]
    S = storage_2d[1:N,1]   # historical daily storage levels [m3]
    q = dailyflow_to_hourly(Q, T); # hourly system net flow (constant over each day) [m3]
    alpha_norm_w = repeat(alpha_norm, N) # replicate available solar capacity for every day of the sim. 
    alpha = alpha_norm_w[1:T*N]

    return L, U, S, q, alpha
end
