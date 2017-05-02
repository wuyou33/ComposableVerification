import numpy as np
import scipy.linalg as linalg
# import timeit
from timeit import Timer
# preprosseseing to generate candidate A's
# t = timeit.Timer("print 'main statement'", "print 'setup'")
def generate_A(n, number_of_instance):
    Q = np.eye(n)
    not_stable = True
    num_A=0
    brutal_time=0
    partitioned_time=0
    while (not_stable or num_A<=number_of_instance-1):
        not_stable = True
        A=np.random.randn(n,n)
        # only interested in the eigenvalues but not the eigenvector
        eig_A = np.linalg.eig(A)[0]
        eig_A.reshape(n,1)
        # print np.real(eig_A)
        # print np.all(np.real(eig_A)<0)
        if np.all(np.real(eig_A)<0):
            not_stable = False
            num_A+=1
            brutal_time+=brutal(A,Q)
            # partitioned_time+=partitioned(A,num_parts)
    print num_A
    print brutal_time
    return brutal_time, partitioned_time


def brutal(A,Q):
    # for a particular A, directly solve 
   t = Timer(lambda: linalg.solve_lyapunov(A,Q))
   single_run_time = t.timeit(number=1)
   print single_run_time
   return single_run_time

def partitioned(A, n, num_parts):
    assert n%num_parts==0
    sub_size = n/num_parts
    A_horizontal=[]
    A_vertical=[]
    for i in len(num_parts):
        for j in enumerate(num_parts):
            A_horizontal.append(A[i*sub_size:(i+1)*sub_size,j*sub_size:(j+1)*sub_size])
        A_vertical.append(A_horizontal)
    print A_horizontal
    print A_vertical

    return A_horizontal, A_vertical
    # X = np.matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))

partitioned(np.eye(60),60,30)
# generate_A(10, 10)