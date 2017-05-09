T = function psd_generator(n)
    B=randn(n,n)
    T=B*B'
end

A = randn(2,2)
A12=randn(2,2)
P = psd_generator(2)
q=psd_generator(2)


# print "eig of P"
println(eig(P)[1])


println(eig(q)[1])

println(eig([P eye(2);q eye(2)])[1])
println(eig(P*q)[1])


# print 'trace of p then q'
# println(trace(P))


# println(trace(q))

# println(trace(P*q))
# # W=randn(1)^2
# # Z=randn(1)^2

# PAAP=P*A*A'*P
# qaaq=A12'*q*q*A12

# print(trace(P*A*q*A12))

# sqrtm_p=sqrtm(PAAP)
# sqrtm_q=sqrtm(qaaq)
# print(trace(sqrtm_p*sqrtm_q))


# show(trace(sqrtm_p*sqrtm_q))


