include("dataload.jl")
include("monte-carlo-functions.jl")

# Load in 2022 - 2023 data & subset for month and year
daily, alpha, RTP = fullsim_dataload();
inflow = daily.inflow
solar_data = alpha[:,1]
price = RTP.MW

# Set data to transform 
data = price
min_val, max_val = minimum(data), maximum(data)

# Simulation parameters
n = length(data)  # Number of time steps
phi = 0.8  # AR(1) coefficient
target_mae = 0.05 # Desired MAE

## If inflow or price
# Normalize data
norm_data = min_max_normalize(data)

# Add noise
simulated_data, ar1_noise = add_noise(n, norm_data, phi, target_mae)

## If inflow or price
# Invert min-max normalization
recovered_data = inverse_minmax(simulated_data, min_val, max_val)

## If solar: Enforce non negative solar capacity
# simulated_data[simulated_data .< 0] .= 0

# Plot results
plot(1:n, recovered_data, label="Simulated Data with AR(1) Errors", color=:orange, linewidth=2)
plot!(1:n, data, label="Original Data", linewidth=2, color=:blue)
# plot!(1:n, ar1_noise, label="AR(1) Noise", linewidth=2, linestyle=:dash)
plot!(legend=:outerbottom)
