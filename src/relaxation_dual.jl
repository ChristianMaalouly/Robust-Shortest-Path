using JuMP
using CPLEX

function relaxationDual(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit=60)

    # JuMP model with the CPLEX optimizer, number of threads, and time limit set
    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 2, "CPX_PARAM_TILIM" => time_limit))
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output

    # Decision variables
    @variable(model, 0 <= x[i=1:n, j=1:n] <= 1)
    @variable(model, 0 <= y[v=1:n] <= 1)
    @variable(model, t1 >= 0)
    @variable(model, t2 >= 0)
    @variable(model, z[i=1:n, j=1:n] >= 0)
    @variable(model, z_prime[v=1:n] >= 0)

    # Objective function
    @objective(model, Min, sum(d[i, j] * x[i, j] for (i, j) in A) + d1 * t1 + sum(D[i, j] * z[i, j] for (i, j) in A))

    # Constraints
    @constraint(model, sum(x[s, i] for i in δplus[s]) - sum(x[i, s] for i in δminus[s]) == 1)
    @constraint(model, sum(x[t, i] for i in δplus[t]) - sum(x[i, t] for i in δminus[t]) == -1)
    @constraint(model, [v in setdiff(1:n, [s, t])], sum(x[i, v] for i in δminus[v]) == sum(x[v, j] for j in δplus[v]))
    @constraint(model, [v in setdiff(1:n, [s, t])], y[v] == sum(x[v, i] for i in δplus[v]))
    @constraint(model, y[s] == 1)
    @constraint(model, y[t] == 1)
    @constraint(model, sum(p[v] * y[v] + 2 * z_prime[v] for v = 1:n) + t2 * d2 <= S)
    @constraint(model, [(i, j) in A], t1 + z[i, j] >= d[i, j] * x[i, j])
    @constraint(model, [v = 1:n], t2 + z_prime[v] >= ph[v] * y[v])


    # Solve the optimization problem
    optimize!(model)

    # Results
    if termination_status(model) == MOI.OPTIMAL
        #println("Optimal solution found.")
        #println("Objective value: ", objective_value(model))
        # Access the values of decision variables if needed
        # x_optimal = value.(x)
        # y_optimal = value.(y)
        # t1_optimal = value(t1)
        # t2_optimal = value(t2)
        # z_optimal = value.(z)
        # z_prime_optimal = value.(z_prime)
        return true, objective_value(model), value.(x)
    else
        #println("No optimal solution found within the time limit.")
        if has_values(model)
            #println("Best solution found:")
            #println("Objective value: ", objective_value(model))
            # Access the values of decision variables if needed
            # x_optimal = value.(x)
            # y_optimal = value.(y)
            # t1_optimal = value(t1)
            # t2_optimal = value(t2)
            # z_optimal = value.(z)
            # z_prime_optimal = value.(z_prime)
            return false, objective_value(model), value.(x)
        else
            return false, Inf, []
        end
    end
end