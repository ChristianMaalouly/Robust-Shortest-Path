using JuMP
using CPLEX
using DataStructures

epsilon = 1e-5

function branchAndCut(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit=60, use_heuristic=false)

    sp1_time = 0
    sp2_time = 0

    # JuMP model with the CPLEX optimizer, number of threads, and time limit set
    model = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_THREADS" => 1, "CPX_PARAM_TILIM" => time_limit))
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

    function my_callback_function(cb_data::CPLEX.CallbackContext, context_id::Clong)
        # If the callback is called because a feasible integer solution is found
        if context_id == CPX_CALLBACKCONTEXT_CANDIDATE
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            x_current = callback_value.(cb_data, x)
            y_current = callback_value.(cb_data, y)
            z_current = callback_value(cb_data, z)

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
            if round(z_current, digits=2) != round(sum(d[i, j] * (1 + delta1_optimal[i, j]) * x_current[i, j] for (i, j) in A), digits=2)
                #println("Subproblem 1 not verified. Adding constraint.")
                #println("Current z value: ", z_current, " Subproblem 1 value: ", sum(d[i, j] * (1 + delta1_optimal[i, j]) * x_current[i, j] for (i, j) in A))
                newConstraint = @build_constraint(z >= sum(d[i, j] * (1 + delta1_optimal[i, j]) * x[i, j] for (i, j) in A))
                MOI.submit(model, MOI.LazyConstraint(cb_data), newConstraint)
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
            if S < sum((p[i] + delta2_optimal[i] * ph[i]) * y_current[i] for i in 1:n) - epsilon
                #println("Subproblem 2 not verified. Adding constraint.")
                #println("Current S value: ", S, " Subproblem 2 value: ", sum((p[i] + delta2_optimal[i] * ph[i]) * y_current[i] for i in 1:n) - epsilon)
                newConstraint = @build_constraint(sum((p[i] + delta2_optimal[i] * ph[i]) * y[i] for i in 1:n) <= S)
                MOI.submit(model, MOI.LazyConstraint(cb_data), newConstraint)
            end
        end

    end

    # Add the callback to the model
    MOI.set(model, CPLEX.CallbackFunction(), my_callback_function)

    # Solve the optimization problem
    optimize!(model)

    # Results
    if termination_status(model) == MOI.OPTIMAL
        #println("Time spent in subproblem 1: ", sp1_time)
        #println("Time spent in subproblem 2: ", sp2_time)
        #println("Optimal solution found.")
        #println("Objective value: ", round(objective_value(model), digits=2))
        # Access the values of decision variables if needed
        #x_optimal = value.(x)
        #y_optimal = value.(y)
        #z_optimal = value.(z)
        return true, objective_value(model), value.(x)
    else
        #println("No optimal solution found within the time limit.")
        if has_values(model)
            #println("Time spent in subproblem 1: ", sp1_time)
            #println("Time spent in subproblem 2: ", sp2_time)
            #println("No optimal solution found within the time limit.")
            #println("Best solution found:")
            #println("Objective value: ", round(objective_value(model), digits=2))
            # Access the values of decision variables if needed
            #x_optimal = value.(x)
            #y_optimal = value.(y)
            #z_optimal = value.(z)
            return false, objective_value(model), value.(x)
        else
            return false, Inf, []
        end
    end

end