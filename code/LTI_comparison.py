import numpy as np
import scipy.linalg as linalg
from timeit import Timer
from cvxpy import *

def lumped(n,blk_size):
    P=Semidef(n,name='lumped_P')
    Q=np.eye(n)
    A=np.random.randn(n,n)
    blks=n/blk_size
    obj = 0 
    final_exp = P*A+A.transpose()*P
    constraints =[
    ]
    prob = cvx.Problem(obj,constraints)

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

lumped(6,3)
# generate_A(10, 10)