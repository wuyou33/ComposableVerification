using JuMP
using Mosek



# generate candiates that admit blkdiag form 
# solutions
function generate_candidate_A(n,num_of_instances)
    Q=eye(n)
    A=randn(n,n)

end

# solves the blkdiag form optimization
# if feasible, returns value, timeit, and stores A

function blkdiag_optimization(A,blk_size)
end

# with the scaling facotors fixed, sovles the low-dimensial optimization, and
# returns the time, the solutions
function LTI_test(partitions=4)
    A = randn(60,60)
    tic()
    P = lyap(A', eye(60))
    toc()
# the composible session
    n1 = size(A,1)
    n2 = size(A,2)
    @assert(n1==n2)
# a method to store the partitions(n)
# test 
end


function composed_test()
    n1 = size(A,1)
    n2 = size(A,2)
    @assert(n1==n2)
end

# LTI_test()
# composed_test()
# adding plotting functionality