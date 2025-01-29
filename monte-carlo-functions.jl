using Random
using Statistics
using Plots

# Function to generate AR(1) Gaussian noise
function generate_ar1_noise(phi, sigma, n)
    noise = zeros(Float64, n)
    # Generate random vector from normal distribution with mean 0 and std 1
    # Scale by tuned sigma
    e = randn(n) * sigma  # White noise
    for t in 2:n
        noise[t] = phi * noise[t-1] + e[t]
    end
    return noise
end

# Function to adjust σ for target MAE
function calibrate_sigma(target_mae, phi, n; tol=1e-3, max_iter=100)
    sigma = target_mae  # Initial guess for σ
    for i in 1:max_iter
        errors = generate_ar1_noise(phi, sigma, n)
        current_mae = mean(abs.(errors))
        if abs(current_mae - target_mae) < tol
            return sigma
        end
        sigma *= target_mae / current_mae  # Adjust σ proportionally
    end
    error("Failed to converge to target MAE")
end

function add_noise(n, data, phi, MAE)
   
    # Calibrate σ
    sigma = calibrate_sigma(target_mae, phi, n)

    # Generate AR(1) noise
    ar1_noise = generate_ar1_noise(phi, sigma, n)

    # Simulate solar forecast data with errors
    simulated_data = data + ar1_noise

    # println("Simulated noise MAE: ", mean(abs.(ar1_noise)))

    return simulated_data, ar1_noise
end
