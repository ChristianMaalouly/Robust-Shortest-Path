using JuMP
using CPLEX
using DataStructures

epsilon = 1e-5

function static(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit=60, use_heuristic=false)

    # JuMP model with the CPLEX optimizer, number of threads, and time limit set
    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 1, "CPX_PARAM_TILIM" => time_limit))
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output

    # Decision variables
    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, y[1:n], Bin)

    # Objective function
    @objective(model, Min, sum(d[i, j] * x[i, j] for (i, j) in A))

    # Constraints
    @constraint(model, sum(x[s, i] for i in δplus[s]) - sum(x[i, s] for i in δminus[s]) == 1)
    @constraint(model, sum(x[t, i] for i in δplus[t]) - sum(x[i, t] for i in δminus[t]) == -1)
    @constraint(model, [v in setdiff(1:n, [s, t])], sum(x[i, v] for i in δminus[v]) == sum(x[v, j] for j in δplus[v]))
    @constraint(model, [v in setdiff(1:n, [s, t])], y[v] == sum(x[v, i] for i in δplus[v]))
    @constraint(model, y[s] == 1)
    @constraint(model, y[t] == 1)
    @constraint(model, sum(p[i] * y[i] for i in 1:n) <= S)

    # Solve the optimization problem
    optimize!(model)

    # Results
    if termination_status(model) == MOI.OPTIMAL
        return true, objective_value(model), value.(x)
    else
        if has_values(model)
            return false, objective_value(model), value.(x)
        else
            return false, Inf, []
        end
    end

end