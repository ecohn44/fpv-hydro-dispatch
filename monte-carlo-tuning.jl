include("monte-carlo-functions.jl")
include("dataload.jl")

inflow = daily.inflow
solar_data = alpha[:,1]
price = RTP.MW

# Set data to transform 
data = solar_data
min_val, max_val = minimum(data), maximum(data)
norm_data = min_max_normalize(data) # For inflow and price
recovered_data = []
simulated_data = []
mapes = [0.05, 0.1, 0.15, 0.2]

# Simulation parameters
n = length(data)  # Number of time steps
phi = 0.8  # AR(1) coefficient
# target_mae = 0.1 # Desired MAE
max_iter=100
tol=1e-3

### Calibrate Sigma for target MAE
for target_mae in mapes
    sigma = target_mae  # Initial guess for σ
    for i in 1:max_iter
        noise = generate_ar1_noise(phi, sigma, n)
        simulated_data = data + noise
        #simulated_data = norm_data + noise
        #recovered_data = inverse_minmax(simulated_data, min_val, max_val)
        
        # current_mae = mean(abs.(errors))
        #current_mae = mean(abs.((recovered_data-data)./(abs.(data).+1)))
        current_mae = mean(abs.((simulated_data-data)./(abs.(data).+.5)))
        
        if abs(current_mae - target_mae) < tol
            break
        end
        sigma *= target_mae / current_mae  # Adjust σ proportionally
    end
    #println("Simulated MAPE: ", mean(abs.((recovered_data-data)./(abs.(data).+1))))
    println("Simulated MAPE: ", mean(abs.((simulated_data-data)./(abs.(data).+.5))))
    println("Sigma: ", sigma)
end

# Add noise
# simulated_data, ar1_noise = add_noise(n, norm_data, phi, target_mae)

## If inflow or price
# Invert min-max normalization
# recovered_data = inverse_minmax(simulated_data, min_val, max_val)

## If solar: Enforce non negative solar capacity
#simulated_data[simulated_data .< 0] .= 0

# Plot results
plot(1:n, simulated_data, label="Simulated Data with AR(1) Errors", color=:orange, linewidth=2)
plot!(1:n, data, label="Original Data", linewidth=2, color=:blue)
#plot!(1:n, noise, label="AR(1) Noise", linewidth=2, linestyle=:dash)
plot!(legend=:outerbottom)
