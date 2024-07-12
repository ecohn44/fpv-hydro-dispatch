import numpy as np

# Define parameters 
T = 10 # number of periods
U = 100 # Cross-period coupling constraint (water contract)
f = lambda x: x**2

# Initialize variables
x = np.zeros(T)
lambda_ = 0.1 # intial lagrangian multiplier
alpha = 0.01 # learning rate

# Lagrangian relaxation iterations
for iter in range(1000):
    # single period optimization 
    for t in range(T):
        x[t] = np.argmin(f(x[t]) + lambda_*x[t])

    # update lagrangian multiplier
    constraint_violation = np.sum(x) - U # water contract constraint in equality formulation 
    lambda_ = lambda_ + alpha*constraint_violation

    # check for convergence 
    if abs(constraint_violation) < 1e-6:
        break

print("Optimal x: ", x)