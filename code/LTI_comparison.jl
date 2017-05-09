using JuMP
using Mosek

# generate candiates that admit blkdiag form 
# solutions
function generate_candidate_A(n,num_of_instances)
    Q=eye(n)
    A=randn(n,n)
    blkdiag_optimization(A,n,blk_size)
    return A, Q
end

# solves the blkdiag form optimization
# if feasible, returns value, timeit, and stores A
function blkdiag_optimization(A,n,blk_size)
    blks=n/blk_size
    solver = MosekSolver()
    m = Model(solver=solver)
    @variable(m,P[1:n,1:n],SDP)
    @constraint(m, constr[i=1:n,j=1:n; j<=(floor((i+blk_size-1)/blk_size)-1)*blk_size||j>(floor((i+blk_size-1)/blk_size))*blk_size], P[i,j] == 0)
    @SDconstraint(m::Model, P*A+A'*P<=-1e-6*eye(n))
    # addinfocallback(m, infocallback, when = :Intermediate)
    solve(m)
    println(MathProgBase.getsolution(m))
end

function infocallback(cb)
    solinfo= MathProgBase.getsolution(cb)
    println(solinfo)
    # push!(bbdata, NodeData(time(),node,obj,bestbound))
end
# addinfocallback(m, infocallback, when = :Intermediate)



function test_setting_zero(n,blk_size)
    P=ones(n,n)
    for i =1:n
        for j = 1:n
            if (j<=(floor((i+blk_size-1)/blk_size)-1)*blk_size|| j>(floor((i+blk_size-1)/blk_size))*blk_size)
                P[i,j] = 0            
            end
        end
    end
    # println(P)
end

function decoupled_LMI(A,blk_size)
    solver = MosekSolver()
    m = Model(solver=solver)
    @variables m begin
    P1[1:blk_size,1:blk_size],SDP
    P2[1:blk_size,1:blk_size],SDP
end
# P=[P1 spzeros(blk_size,blk_size);spzeros(blk_size,blk_size) P2]
@SDconstraint(m::Model, P*A+A'*P<=-eye(60))
solve(m)
end

function composed_test()
    n1 = size(A,1)
    n2 = size(A,2)
    @assert(n1==n2)
end

# test_setting_zero(8,4)
t=@time(blkdiag_optimization(-eye(60,60),60,3))
# @time(generate_candidate_A(3,100))

# adding plotting functionality