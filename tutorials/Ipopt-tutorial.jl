using JuMP, Ipopt

m = Model(Ipopt.Optimizer)

@variable(m, 0 <= p,  base_name="Quantities of pizzas")
@variable(m, 0 <= s,  base_name="Quantities of sandwiches")

@constraint(m, budget,     10p + 4s <=  80 )

@objective(m, Max, 100*p - 2*p^2 + 70*s - 2*s^2 - 3*p*s)

optimize!(m)

status = termination_status(m)

if (status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT) && has_values(m)
    if (status == MOI.OPTIMAL)
        println("** Problem solved correctly **")
    else
        println("** Problem returned a (possibly suboptimal) solution **")
    end
    println("- Objective value : ", objective_value(m))
    println("- Optimal solutions:")
    println("pizzas: $(value.(p))")
    println("sandwitches: $(value.(s))")
    println("- Dual (budget): $(dual.(budget))")
else
    println("The model was not solved correctly.")
    println(status)
end