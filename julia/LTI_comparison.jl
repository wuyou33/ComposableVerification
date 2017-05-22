using JuMP
using Mosek

# generate candiates that admit blkdiag form 
# solutions
function generate_candidate_A(n,num_of_instances)
    Q=eye(n)
    A=randn(n,n)
    lumped_optimization(A,n,blk_size)
    return A, Q
end

function write_data(A,brutal_time,opt_scaling,fixed_scaling)
    open("data.csv","w") do fp
    println(fp, "A, brutal_time, opt_scaling, fixed_scaling")
    for bb in bbdata
        println(fp, A, ",", brutal_time, ",",
            opt_scaling, ",", fixed_scaling)
    end
end
end

# solves the blkdiag form optimization
# if feasible, returns value, timeit, and stores A
function lumped_optimization(A,n,blk_size)
    Q=eye(n)
    A=randn(n,n)

    blks=n/blk_size
    m = Model(solver=MosekSolver(MSK_IPAR_LOG=0))
    @variable(m,P[1:n,1:n],SDP)
    @constraint(m, constr[i=1:n,j=1:n; j<=(floor((i+blk_size-1)/blk_size)-1)*blk_size||j>(floor((i+blk_size-1)/blk_size))*blk_size], P[i,j] == 0)
    @SDconstraint(m::Model, P*A+A'*P<=-1e-6*eye(n))
    # addinfocallback(m, infocallback, when = :Intermediate)
    solve_status=solve(m)
    if solve_status==:Optimal
        t = getsolvetime(m)
        scaled_time = scaled_decouple(A,n,blk_size)
        fixed_time = fixed_decouple(A,n,blk_size)
        # println(t)
    end
end



function lumped_optimization_2(A,n,blk_size)
    Q=eye(n)
    A=randn(n,n)

    blks=n/blk_size
    m = Model(solver=MosekSolver(MSK_IPAR_LOG=0))
    @variable(m,P[1:blk_size,1:n])
    # @constraint(m, constr[i=1:n,j=1:n; j<=(floor((i+blk_size-1)/blk_size)-1)*blk_size||j>(floor((i+blk_size-1)/blk_size))*blk_size], P[i,j] == 0)
    for i =1:blks
        @SDconstraint(m::Model, P[:,i:i+1]>=1e-6*eye(1)) 
    end

    @SDconstraint(m::Model, P*A+A'*P<=-1e-6*eye(n))
    # addinfocallback(m, infocallback, when = :Intermediate)
    solve_status=solve(m)
    # if solve_status==:Optimal
    #     t = getsolvetime(m)
    #     scaled_time = scaled_decouple(A,n,blk_size)
    #     fixed_time = fixed_decouple(A,n,blk_size)
    #     # println(t)
    # end
end



function scaled_decouple(A,n,blk_size)
    blks=n/blk_size
    m = Model(solver=MosekSolver(MSK_IPAR_LOG=0))
    @variable(m,P[1:blk_size,1:n])
    @variable(m,M[1:blk_size,1:blk_size],SDP)

    for i = 1:blks
        @SDconstraint(m::Model, P[:,(i-1)*blk_size+1:i*blk_size]>=1e-6*eye(n)) 
    end
end

function fixed_decouple(A,n,blk_size)
end

function test_setting_sparsity(n,blk_size)
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
lumped_optimization_2(-eye(6,6),6,3)
# Riccati_based(-eye(6,6),6,3)

# @time(generate_candidate_A(3,100))


# adding plotting functionality


# utility functions

function infocallback(cb)
    solinfo= MathProgBase.getsolution(cb)
    println(solinfo)
    # push!(bbdata, NodeData(time(),node,obj,bestbound))
end
# addinfocallback(m, infocallback, when = :Intermediate)



