using JuMP
using CPLEX
using DataStructures

epsilon = 1e-5

function enhancedCuttingPlanes(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit=60, use_heuristic=false)

    sp1_time = 0
    sp2_time = 0
    start_time = time()
    optimal_solution = n * maximum(d) * maximum(D)
    optimal_x = []

    # JuMP model with the CPLEX optimizer, number of threads, and time limit set
    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 2, "CPX_PARAM_TILIM" => time_limit))
    set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output

    # Decision variables
    @variable(model, x[1:n, 1:n], Bin)
    @variable(model, y[1:n], Bin)
    @variable(model, z)

    # Objective function
    @objective(model, Min, z)

    # Constraints
    @constraint(model, z >= sum(d[i, j] * x[i, j] for (i, j) in A))

    @constraint(model, sum(x[s, i] for i in δplus[s]) - sum(x[i, s] for i in δminus[s]) == 1)
    @constraint(model, sum(x[t, i] for i in δplus[t]) - sum(x[i, t] for i in δminus[t]) == -1)
    @constraint(model, [v in setdiff(1:n, [s, t])], sum(x[i, v] for i in δminus[v]) == sum(x[v, j] for j in δplus[v]))
    @constraint(model, [v in setdiff(1:n, [s, t])], y[v] == sum(x[v, i] for i in δplus[v]))
    @constraint(model, y[s] == 1)
    @constraint(model, y[t] == 1)

    @constraint(model, sum(p[i] * y[i] for i in 1:n) <= S)


    updated = true
    # Solve MP, then check SP1 and SP2 and repeat if needed
    while updated && time() - start_time < time_limit
        updated = false
        sp1 = false

        # Solve MP
        optimize!(model)

        if !has_values(model)
            break
        end

        # Get solution values
        x_current = value.(x)
        y_current = value.(y)
        z_current = value(z)

        #println("Current z value: ", z_current, " Optimal solution: ", optimal_solution)
        if z_current > optimal_solution
            break
        end

        # SP1
        delta1_optimal = zeros(n, n)
        start_time_sp1 = time()

        if use_heuristic
            # Heuristic
            # Create a priority queue with all distances used in the solution
            distances = PriorityQueue{Tuple{Int,Int},Float64}(Base.Order.Reverse)
            for (i, j) in A
                if x_current[i, j] > 0
                    enqueue!(distances, (i, j) => d[i, j])
                end
            end

            # Increase the distances as much as possible based on the limit d1
            d1_used = 0

            while d1_used < d1 && !isempty(distances)
                i, j = dequeue!(distances)
                if d1_used + D[i, j] <= d1
                    d1_used += D[i, j]
                    delta1_optimal[i, j] = D[i, j]
                else
                    delta1_optimal[i, j] = d1 - d1_used
                    d1_used = d1
                end
            end
        else
            # SP1 model
            model_sp1 = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 1))
            set_optimizer_attribute(model_sp1, "CPX_PARAM_SCRIND", 0) # Remove the solver output

            @variable(model_sp1, delta1_sp[1:n, 1:n] >= 0)
            @objective(model_sp1, Max, sum(d[i, j] * (1 + delta1_sp[i, j]) * x_current[i, j] for (i, j) in A))
            @constraint(model_sp1, sum(delta1_sp[i, j] for (i, j) in A) <= d1)
            @constraint(model_sp1, delta1_sp .<= D)

            optimize!(model_sp1)
            delta1_optimal = value.(delta1_sp)
        end
        # Sum of time taken by all SP1 solutions
        sp1_time += time() - start_time_sp1

        # Verify if truly optimal, otherwise add constraint
        new_opt_sol = round(sum(d[i, j] * (1 + delta1_optimal[i, j]) * x_current[i, j] for (i, j) in A), digits=2)
        if round(z_current, digits=2) != new_opt_sol
            """
            # Print path and uncertainty result to analyse
            println("Current z value: ", z_current, " Subproblem 1 value: ", new_opt_sol)
            for (i, j) in A
                if x_current[i, j] > 0
                    println(i, " ", j, " ", D[i, j], " ", d[i, j], " ", x_current[i, j], " ", delta1_optimal[i, j])
                end
            end
            println("--------")
            println(sum(d[i, j] * (1 + D[i, j]) for (i, j) in A if x_current[i, j] > 0.5))
            """
            #println("Subproblem 1 not verified. Adding constraint.")
            @constraint(model, z >= sum(d[i, j] * (1 + delta1_optimal[i, j]) * x[i, j] for (i, j) in A))
            sp1 = true
            updated = true
        end

        # SP2
        start_time_sp2 = time()
        delta2_optimal = zeros(n)

        if use_heuristic
            # Heuristic
            # Create a priority queue with all weights used in the solution
            weights = PriorityQueue{Int,Float64}(Base.Order.Reverse)
            for i in 1:n
                if y_current[i] > 0
                    enqueue!(weights, i => ph[i])
                end
            end

            # Increase the distances as much as possible based on the limit d1
            full_iterations = min(d2 ÷ 2, length(weights))
            remainder = d2 % 2
            for _ in 1:full_iterations
                i = dequeue!(weights)
                delta2_optimal[i] = 2
            end
            if remainder > 0 && !isempty(weights)
                i = dequeue!(weights)
                delta2_optimal[i] = remainder
            end

        else
            # SP2 model
            model_sp2 = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 1))
            set_optimizer_attribute(model_sp2, "CPX_PARAM_SCRIND", 0) # Remove the solver output

            @variable(model_sp2, 0 <= delta2_sp[1:n] <= 2)
            @objective(model_sp2, Max, sum((p[i] + delta2_sp[i] * ph[i]) * y_current[i] for i in 1:n))
            @constraint(model_sp2, sum(delta2_sp) <= d2)

            optimize!(model_sp2)
            delta2_optimal = value.(delta2_sp)
        end
        # Sum of time taken by all SP2 solutions
        sp2_time += time() - start_time_sp2

        # Verify if feasible to update optimal solution, otherwise add constraint
        #println("optimal solution: ", optimal_solution, " new optimal solution: ", new_opt_sol)
        if S < sum((p[i] + delta2_optimal[i] * ph[i]) * y_current[i] for i in 1:n) + epsilon
            #println("Subproblem 2 not verified. Adding constraint.")
            #println("Current S value: ", S, " Subproblem 2 value: ", sum((p[i] + delta2_optimal[i] * ph[i]) * y_current[i] for i in 1:n) - epsilon)
            @constraint(model, sum((p[i] + delta2_optimal[i] * ph[i]) * y[i] for i in 1:n) - epsilon <= S)
            updated = true
        else
            if sp1
                @constraint(model, sum(d[i, j] * (1 + D[i, j]) * x[i, j] for (i, j) in A) <= sum(d[i, j] * (1 + D[i, j]) for (i, j) in A if x_current[i, j] > 0.5) - 1)
            end
            if new_opt_sol < optimal_solution
                optimal_solution = new_opt_sol
                optimal_x = x_current
            end
        end

        """
        if updated && time() - start_time < time_limit
            set_start_value.(x, x_current)
            set_start_value.(y, y_current)
            set_start_value(z, z_current)
        end
        """
    end

    """
    # Print the path found to analyse results
    println("final result")
    x_current = optimal_x
    for (i, j) in A
        if x_current[i, j] > 0
            println(i, " ", j, " ", D[i, j], " ", d[i, j])
        end
    end
    """

    # Results
    if time() - start_time < time_limit
        #println("Time spent in subproblem 1: ", sp1_time)
        #println("Time spent in subproblem 2: ", sp2_time)
        #println("Optimal solution found.")
        #println("Objective value: ", optimal_solution)
        # Access the values of decision variables if needed
        #x_optimal = optimal_x

        return true, optimal_solution, optimal_x
    else
        #println("Time spent in subproblem 1: ", sp1_time)
        #println("Time spent in subproblem 2: ", sp2_time)
        #println("No optimal solution found within the time limit.")
        #println("Best solution found:")
        #println("Objective value: ", optimal_solution)
        # Access the values of decision variables if needed
        #x_optimal = optimal_x

        return false, optimal_solution, optimal_x
    end

end