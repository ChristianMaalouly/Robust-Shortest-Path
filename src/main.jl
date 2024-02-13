include("read_data.jl")
include("dual.jl")
include("cutting_planes.jl")
include("branch_and_cut.jl")
include("relaxation_dual.jl")
include("enhanced_cutting_planes.jl")
include("static.jl")

# Heuristic method
function heuristic(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit)

    # Prpare the new set of edges
    new_A = []

    # Start timer
    start_time = time()

    # Solved the relaxed problem
    optimal, solution, x = relaxationDual(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit)
    if optimal

        # Add the edges to the new set
        # Update the distances in a temporary list for more relaxed solutions
        new_d = deepcopy(d)
        for i in 1:n
            for j in 1:n
                if x[i, j] > 0
                    push!(new_A, (i, j))
                    new_d[i, j] = 2 * d[i, j]
                    #println(i, " ", j, " ", x[i, j], " ", d[i, j], " ", D[i, j])
                end
            end
        end

        # Update the time passed
        time_passed = time() - start_time
        if time_passed < time_limit

            # Solve the relaxed problem again
            optimal, solution, x = relaxationDual(n, s, t, δplus, δminus, A, new_d, d1, d2, D, p, ph, S, time_limit - time_passed)
            if optimal
                # Add the edges to the new set
                for i in 1:n
                    for j in 1:n
                        if x[i, j] > 0
                            if (i, j) ∉ new_A
                                push!(new_A, (i, j))
                                #println(i, " ", j, " ", x[i, j], " ", d[i, j], " ", D[i, j])
                            end
                        end
                    end
                end

                # Get the neighbor sets for the new set of edges
                new_δplus, new_δminus = getNeighbours(n, new_A)

                # Update the time passed
                time_passed = time() - start_time
                if time_passed < 160

                    # Solve with the dualisation method using the new edges with the original distances
                    return dualisationMethod(n, s, t, new_δplus, new_δminus, new_A, d, d1, d2, D, p, ph, S, time_limit - time_passed)
                end
            else
                return false, Inf, []
            end
        end
    else
        return false, Inf, []
    end
end

methods = ["Dualisation", "Cutting planes", "Cutting planes with heuristic", "Enhanced cutting planes", "Branch and cut", "Branch and cut with heuristic", "Heuristic", "Relaxed dualisation", "Static problem"]
println("main() function takes 3 integers, \nfirst for the method, second for the instance, and third for the time limit in seconds\n")
println("for the method: \n1: dualisation\n2: cutting planes\n3: cutting planes with heuristic\n4: enhanced cutting planes\n5: branch and cut\n6: branch and cut with heuristic\n7: heuristic\n8: relaxed dualisation\n9: static problem\n")
println("for the instances: \n1 to 41: BAY instanced of size 20 to 2500\n42 to 82: COL instanced of size 20 to 2500\n83 to 123: NY instanced of size 20 to 2500\n0: for all instances")
function main(method, instance, time_limit)
    n, s, t, S, d1, d2, p, ph, A, d, D = read_data(dataFiles[instance])
    δplus, δminus = getNeighbours(n, A)
    println("Instance: ", dataFiles[instance])


    println(methods[method])
    start_time = time()
    if method == 1
        optimal, solution, x = dualisationMethod(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit)
    elseif method == 2
        optimal, solution, x = cuttingPlanes(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit)
    elseif method == 3
        optimal, solution, x = cuttingPlanes(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit, true)
    elseif method == 4
        optimal, solution, x = enhancedCuttingPlanes(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit, true)
    elseif method == 5
        optimal, solution, x = branchAndCut(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit)
    elseif method == 6
        optimal, solution, x = branchAndCut(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit, true)
    elseif method == 7
        optimal, solution, x = heuristic(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit)
    elseif method == 8
        optimal, solution, x = relaxationDual(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit)
    elseif method == 9
        optimal, solution, x = static(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, time_limit)
    else
        println("bad input for method")
        return
    end

    if optimal
        println("Optimal solution found!")
    else
        println("No optimal solution found within the time limit!")
    end
    println("Solution value: ", solution)
    println("Time: ", time() - start_time)

    println("")
end



#for i in 1:length(dataFiles)
#    result_file = open("results_cut_heuristic.txt", "a")
#    println(result_file, "Running for ", dataFiles[i])
#    n, s, t, S, d1, d2, p, ph, A, d, D = read_data(dataFiles[i])
#    δplus, δminus = getNeighbours(n, A)
#
#    start_time = time()
#    optimal, solution, x = cuttingPlanes(n, s, t, δplus, δminus, A, d, d1, d2, D, p, ph, S, 30, true)
#    if optimal
#        println(result_file, "Optimal solution found!")
#    else
#        println(result_file, "No optimal solution found within the time limit!")
#    end
#    println(result_file, "Solution value: ", solution)
#    #println(result_file, "Solution: ", x)
#    println(result_file, "Time: ", time() - start_time)
#
#    println(result_file, "")
#    println(result_file, "==============================")
#    println(result_file, "")
#    close(result_file)
#end
