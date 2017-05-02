import numpy as np
import scipy.linalg as linalg
import timeit
# preprosseseing to generate candidate A's
t = timeit.Timer("print 'main statement'", "print 'setup'")
def generate_A(n, number_of_instance):
    Q = np.eye(n)
    not_stable = True
    num_A=0
    brutal_time=0
    partitioned_time=0
    while (not_stable or num_A<=number_of_instance-1): 
        A=np.random.randn(n,n)
        # only interested in the eigenvalues but not the eigenvector
        eig_A = np.linalg.eig(A)[0]
        eig_A.reshape(n,1)
        # print np.real(eig_A)
        # print np.all(np.real(eig_A)<0)
        if np.all(np.real(eig_A)<0):
            not_stable = False
            num_A+=1
            brutal(A,Q)
            # brutal_time+=brutal(A,Q)
            # print brutal_time
    return brutal_time, partitioned_time


def brutal(A,Q):
    # for a particular A, directly solve 
   t.timeit(linalg.solve_lyapunov(A,Q))
    # print timeit.timeit('char in text', setup='text = "sample string"; char = "g"')
    # return time_spent

def partitioned(A, num_parts):
    #     X = np.matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))
    return time_spent 

generate_A(6, 10)