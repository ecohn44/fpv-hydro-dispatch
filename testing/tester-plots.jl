using CSV
using DataFrames
using Plots

# Load the two CSV files
relaxed_revenue = CSV.read("testing/relaxed_revenue_gurobi_optwc.csv", DataFrame)
baseline_revenue = CSV.read("testing/baseline_revenue_gurobi_wc.csv", DataFrame)

# Extract the columns (assuming each file has one column)
relaxed_values = relaxed_revenue[:, 1]  # First column
baseline_values = baseline_revenue[:, 1]  # First column

# Create a plot and overlay the two datasets
plot(relaxed_values/1e6, label="Relaxed Revenue", linewidth=2, color=:blue)
plot!(baseline_values/1e6, label="Baseline Revenue", linewidth=2, color=:red)
title!("Comparing Revenue over N with Water Contract")
xlabel!("Horizon Length (day)")
ylabel!("Revenue (\$ Million)")