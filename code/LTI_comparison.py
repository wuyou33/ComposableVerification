import numpy as np
import scipy.linalg as linalg
from timeit import Timer
from cvxpy improt *

def lumped(n,blk_size)
    Q=np.eye(n)
    A=np.randn(n,n)
    blks=n/blk_size
    m = Model(solver=MosekSolver(MSK_IPAR_LOG=0))
    @variable(m,P[1:blk_size,1:n])
        # addinfocallback(m, infocallback, when = :Intermediate)
    solve_status=solve(m)


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