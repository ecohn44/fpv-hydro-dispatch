using Statistics
include("monte-carlo-functions.jl")
include("dataload.jl")

daily, alpha, RTP = fullsim_dataload();
inflow = daily.inflow
# Compute the average profile over columns
solar_data = mean.(eachrow(alpha))
price = RTP.MW

# Set data to transform 
data = solar_data
solar = true

# Apply normalization (redundant for solar)
min_val, max_val = minimum(data), maximum(data)
norm_data = min_max_normalize(data) # For inflow and price
recovered_data = []
simulated_data = []
mapes = [0 , 0.05, 0.1, 0.15, 0.2]

# Simulation parameters
n = length(data)  # Number of time steps
phi = 0.8  # AR(1) coefficient
max_iter=100
tol=1e-3

### Calibrate Sigma for target MAE
for target_mae in mapes
    sigma = target_mae  # Initial guess for σ
    for i in 1:max_iter
        
        noise = generate_ar1_noise(phi, sigma, n)

        if solar 
            simulated_data = data + noise
            current_mae = mean(abs.((simulated_data-data)./(abs.(data).+.5)))
        else 
            simulated_data = norm_data + noise
            recovered_data = inverse_minmax(simulated_data, min_val, max_val)
            current_mae = mean(abs.((recovered_data-data)./(abs.(data).+1)))
        end 
        
        if abs(current_mae - target_mae) < tol
            break
        end
        sigma *= target_mae / current_mae  # Adjust σ proportionally
    end

    if solar 
        println("Simulated MAPE: ", mean(abs.((simulated_data-data)./(abs.(data).+.5))))
        simulated_data[simulated_data .< 0] .= 0
    else 
        println("Simulated MAPE: ", mean(abs.((recovered_data-data)./(abs.(data).+1))))
    end 
    println("Sigma: ", sigma)
end


# Plot results
if solar 
    plot(1:n, simulated_data, label="Simulated Data with AR(1) Errors", color=:orange, linewidth=2)
else
    plot(1:n, recovered_data, label="Simulated Data with AR(1) Errors", color=:orange, linewidth=2)
end 
plot!(1:n, data, label="Original Data", linewidth=2, color=:blue)
#plot!(1:n, noise, label="AR(1) Noise", linewidth=2, linestyle=:dash)
plot!(legend=:outerbottom)
