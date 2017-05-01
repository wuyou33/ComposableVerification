using JuMP

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

LTI_test()
# composed_test()
