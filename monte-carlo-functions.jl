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

function add_noise(sigma, phi, n, data)
   
    # Calibrate σ
    # sigma = calibrate_sigma(target_mae, phi, n)

    # Generate AR(1) noise
    ar1_noise = generate_ar1_noise(phi, sigma, n)

    # Simulate solar forecast data with errors
    simulated_data = data + ar1_noise
    
    # Enforce positivity for solar input
    simulated_data[simulated_data .< 0] .= 0

    println("Simulated MAPE: ", mean(abs.((simulated_data-data)./(abs.(data).+.5))))

    return simulated_data
end

function add_noise_norm(sigma, phi, n, data)

    min_val, max_val = minimum(data), maximum(data)
    norm_data = min_max_normalize(data)
    noise = generate_ar1_noise(phi, sigma, n)
    simulated_data = norm_data + noise
    recovered_data = inverse_minmax(simulated_data, min_val, max_val)
        
    println("Simulated MAPE: ", mean(abs.((recovered_data-data)./(abs.(data).+1))))

    return recovered_data
end
