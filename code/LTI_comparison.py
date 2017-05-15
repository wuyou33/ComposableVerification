import sys
import mosek
from   mosek.fusion import *
import numpy as np

def lumped(n,blk_size,feasibility_only=True, quite=True):

    Q=np.eye(n)
    A=np.random.randn(n,n)
    # A=-np.rand(n)
    assert n%blk_size==0
    blks=n/blk_size
    with Model("lumped") as M:    
        # Setting up the variables
        P = M.variable("P", Domain.inPSDCone(n))

        # structure constraint
        for i in range(n):
            for j in range(n):
                    if j >=((i/blk_size)+1)*blk_size or j< (i/blk_size)*blk_size:
                        M.constraint(P.pick([i],[j]), Domain.equalsTo(0.0))

        # lyapunov constraints
        M.constraint(Expr.neg(Expr.add(Expr.mul(P,A),Expr.mul(A.transpose(),P))),Domain.inPSDCone())


        if not feasibility_only:
            M.objective(ObjectiveSense.Minimize, Expr.dot(P,Q))

        if not quite:
            M.setLogHandler(sys.stdout)

        M.solve()
        status = M.getSolverDoubleInfo("intpntOptStatus")
        # print("solver status " + str(status))
        if status == 1:
            lumped_time = M.getSolverDoubleInfo("intpntTime")

            # print("solver lumped_time " + str(lumped_time))
            # print(P.level())
        else:
            lumped_time = 'inf'    
        return A, lumped_time, P


# lumped(6,4)

# form n codependent LMIs of the origianl problem size to search for good scaling factors
def direct_scaling(A, n, blk_size):
    blks=n/blk_size
    P=[]
    N=[]
    with Model("lumped") as M:    
        # Setting up the variables
        for i in range(blks):
            P.append(M.variable(Domain.inPSDCone(blk_size)))
            N.append(M.variable(Domain.inPSDCone(blk_size)))


# w/ fixed scaling factor 
def fixed_scaling(A,blk_size):
    pass
    

def ricaati_based(A,Q):
    pass
    # t = Timer(lambda: linalg.solve_lyapunov(A,Q))
        # X = np.matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))


def read_data(file, key):
    pass

def write_data(key, value):
    pass





