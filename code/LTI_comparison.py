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

# form n codependent LMIs of the origianl problem size to search for good scaling factors
def direct_scaling(A, n, blk_size, feasibility_only=True, quite=True):
    blks=n/blk_size
    P=[]
    N=[]
    PAAP=[]
    with Model("direct_scaling") as M:    
        # Setting up the variables
        for i in range(blks):
            P.append(M.variable(Domain.inPSDCone(blk_size)))
            N.append(M.variable(Domain.inPSDCone(blk_size)))
        # for i in range(blks):
            # extracting the corresponding A submatrices
            # print(P[i])
            row_slice = Expr.mul(P[i],A[(i)*blk_size:(i+1)*blk_size,:])

            colum_slice = row_slice.transpose()
            # print(row_slice.size())
            if i == 0:
                rowstaced=Expr.vstack(row_slice,Expr.reshape(Expr.zeros((blks-1)*blk_size*n),(blks-1)*blk_size,n))
                columnstaced=Expr.hstack(colum_slice,Expr.reshape(Expr.zeros((blks-1)*blk_size*n),n,(blks-1)*blk_size))
            elif i==n:
                rowstaced=Expr.vstack(Expr.reshape(Expr.zeros((blks-1)*blk_size*n),(blks-1)*blk_size,n),row_slice)
                columnstaced=Expr.hstack(Expr.reshape(Expr.zeros((blks-1)*blk_size*n),n,(blks-1)*blk_size),colum_slice)
            else:
                ontop=Expr.reshape(Expr.zeros((i)*blk_size*n),(i)*blk_size,n)
                onbuttom=Expr.reshape(Expr.zeros((blks-i-1)*blk_size*n),(blks-i-1)*blk_size,n)
                rowstaced=Expr.vstack(ontop,row_slice,onbuttom)

                onLeft=Expr.reshape(Expr.zeros((i)*blk_size*n),n,(i)*blk_size)
                onRight=Expr.reshape(Expr.zeros((blks-i-1)*blk_size*n),n,(blks-i-1)*blk_size)
                columnstaced=Expr.hstack(onLeft,colum_slice,onRight)
            PAAP.append(Expr.add(rowstaced,columnstaced))
        # print(PAAP[0])
        # setting up the constraints, essentially solving for n independent LMIs
        # M.constraint(Expr.neg(member_LMI),Domain.inPSDCone())
direct_scaling(np.random.randn(6,6),6,2)

def construct_blk_diag(N, blk_size):
    zero_blks = Expr.reshape(Expr.zeros(blk_size*blk_size),blk_size,blk_size)

    for i in range(len(N)):
        Expr.hstack(N[i],zero_blks)



    return blk_diag_form 

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





# testing_stuff()