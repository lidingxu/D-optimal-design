# julia jumpsolve.jl instance time_limit solver_name
using JuMP
import Gurobi
import CPLEX
import SCIP
using LinearAlgebra
using MathOptInterface
const MOI = MathOptInterface


struct Problem
    numvars::Int
    dim::Int
    cap::Int
    data::Matrix{Float64}
    epsilon::Float64

    function Problem(numvars::Int, dim::Int, cap::Int, data::Matrix{Float64}, epsilon::Float64)
        problem = new(numvars, dim, cap, data, epsilon)
        return problem
    end

end

# solve:
function solve(problem::Problem, time_limit::Float64, solver::String = "CPLEX")
    epsilon = problem.epsilon
    numvars = problem.numvars
    dim = problem.dim
    cap = problem.cap
    data = problem.data
    print("numvars:", numvars, " dim:", dim, " epsilon", epsilon, "rank ", rank(data), "\n")
    bininds = 1:numvars
    binindsplus = 1:(numvars+1)
    diminds = 1:dim
    if solver == "CPLEX"
        m = Model(CPLEX.Optimizer) 
        set_optimizer_attribute(m, "CPXPARAM_TimeLimit", time_limit) 
    elseif solver == "GUROBI"
        m = Model(Gurobi.Optimizer)
        set_optimizer_attribute(m, "TimeLimit",  time_limit)
    elseif solver == "SCIP"
        m = Model(SCIP.Optimizer)
        set_optimizer_attribute(m, "limits/time",  time_limit)
    else
        return 0
    end



    @variables(m, begin
        x[i in bininds], Bin # 
        Z[i in bininds, j in diminds]
        J[j1 in diminds, j2 in diminds]
        t[i in binindsplus, j in diminds] >= 0
        diagJ[j in diminds]
        obj >= 0
        epsZ[j1 in diminds, j2 in diminds]
    end
    )


    
    # bounds
    xvals = [-0.0, -0.0, 1.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, 1.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, 1.0, -0.0, 1.0]
    @constraints(m, begin
        #[i in bininds], x[i] == xvals[i]
        [j1 in diminds, j2 in (j1+1):dim], J[j1, j2] == 0 # J[j1][j2] = 0, j2 >= j1
        [j in diminds], J[j, j] >= 0 # J[j][j] >= 0
        [j1 in diminds, j2 in j1:dim], sum(data[i, j1]*Z[i, j2] for i in bininds) +  epsilon * epsZ[j1, j2] == J[j1, j2] # sum(AiZi) = j
        [i in bininds, j in diminds], Z[i,j]* Z[i,j] <= 2 * t[i,j] * x[i] # [t[i, j]; x[i]; Z[i, j]] in RotatedSecondOrderCone() # w^2 <= z[2] * s
        [j in diminds],  sum(epsZ[i, j] * epsZ[i, j] for i in diminds)  <= 2 * t[numvars + 1,j]
        [j in diminds],  sum(t[i, j] for i in binindsplus) * 2 <= J[j, j] # sum(tij) <= Jjj
        [j in diminds], diagJ[j] == J[j,j]
        [obj; diagJ] in MOI.GeometricMeanCone(dim + 1) 
        sum(x[i] for i in bininds) <= cap
    end
    )    
    #
    
    @objective(m, Max,  obj)

    #write_to_file(m, "my_file.mps")
    #print(Omega) 
    #print(m)
    #return 0
    
 

    optimize!(m) 

    print([value(x[i]) for i in bininds])

    log_det_from = false
    if log_det_from
        bd = objective_bound(m)
        bd = log(bd) -  2 * log(epsilon)
        sol_val = log(objective_value(m)) -  2 * log(epsilon)
        val =  log(det(transpose(data) * Diagonal([value(x[i]) for i in bininds]) * data + epsilon^2 * I))
        val =  val/dim - 2 * log(epsilon)
        print(val, " ", sol_val, " ", bd, "\n")
    else
        sol_val = objective_value(m) -  epsilon^2
        val =  det(transpose(data) * Diagonal([value(x[i]) for i in bininds]) * data + epsilon^2 * I)^(1.0/dim)
        val =  val - epsilon^2
        print(val, " ", sol_val, "\n")
    end

    return sol_val
end


function readProblem(absolute_path)
    lines = readlines(absolute_path)
    line1 = split.(lines[1], " ")
    numvars = parse(Int, String(line1[1]))
    dim = parse(Int, String(line1[2]))
    data = zeros(Float64, (numvars, dim))
    cap = dim + 2

    i = 1
    for line in lines[2:end]
        line = split.(line, " ")
        j = 1
        for entry in line
            if entry == ""
                break
            end
            data[i, j]  = parse(Float64, String(entry))
            j += 1
        end
        i += 1
    end

    epsilon = 1e-3

    problem = Problem(numvars, dim, cap, data, epsilon)
    return problem
end

function main(args)
    print("welcome\n")
    instance = args[1]
    time_limit =  parse(Float64, args[2])
    solver =   args[3]
    problem = readProblem(instance)
    print("problem loaded\n")
    solve(problem, time_limit, solver)
end

main(ARGS)
