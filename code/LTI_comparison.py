import sys
import mosek
from   mosek.fusion import *
import numpy as np
def construct_blk_diag(ingredients, num_blks, blk_size, index, negate=True):
    blk_diag_form=Expr.reshape(Expr.zeros(0),0,num_blks*blk_size)
    for i in range(num_blks):
        onLeft=Expr.reshape(Expr.zeros((i)*blk_size*blk_size),blk_size,(i)*blk_size)
        onRight=Expr.reshape(Expr.zeros((num_blks-i-1)*blk_size*blk_size),blk_size,(num_blks-i-1)*blk_size)
        if negate and index==i:
            row_ingredient=Expr.hstack(onLeft,Expr.neg(ingredients[i]),onRight)
        else:
            row_ingredient=Expr.hstack(onLeft,ingredients[i],onRight)
        blk_diag_form=Expr.vstack(blk_diag_form,row_ingredient)
    return blk_diag_form


def lumped(n,blk_size,feasibility_only=True, quite=True):
    Q=np.eye(n)
    A=np.random.randn(n,n)
    # A=-np.rand(n)
    assert n%blk_size==0
    num_blks=n/blk_size
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
        else:
            lumped_time = 'inf'    
        return A, lumped_time, P

# form n codependent LMIs of the origianl problem size to search for good scaling factors
def direct_scaling(A, n, blk_size, feasibility_only=True, quite=True):
    num_blks=n/blk_size
    P=[]
    N=[]
    PAAP=[]
    with Model("direct_scaling") as M:    
        # Setting up the variables
        for i in range(num_blks):
            P.append(M.variable(Domain.inPSDCone(blk_size)))
            N.append(M.variable(Domain.inPSDCone(blk_size)))

        # Setting up the constraints, i is the constraint index
        for i in range(num_blks):
            row_slice = Expr.neg(Expr.mul(P[i],A[(i)*blk_size:(i+1)*blk_size,:]))
            colum_slice = row_slice.transpose()

            ontop=Expr.reshape(Expr.zeros((i)*blk_size*n),(i)*blk_size,n)
            onbuttom=Expr.reshape(Expr.zeros((num_blks-i-1)*blk_size*n),(num_blks-i-1)*blk_size,n)
            rowstaced=Expr.vstack(ontop,row_slice,onbuttom)

            onLeft=Expr.reshape(Expr.zeros((i)*blk_size*n),n,(i)*blk_size)
            onRight=Expr.reshape(Expr.zeros((num_blks-i-1)*blk_size*n),n,(num_blks-i-1)*blk_size)
            columnstaced=Expr.hstack(onLeft,colum_slice,onRight)
            row_and_column = Expr.add(rowstaced,columnstaced)
            scaling=construct_blk_diag(N, num_blks, blk_size, i, negate=True)
            PAAP.append(Expr.add(row_and_column,scaling))
        for i in range(num_blks):
            print("add constraint")
        # setting up the constraints, essentially solving for n independent LMIs
            M.constraint(PAAP[i],Domain.inPSDCone())

        if not feasibility_only:
            M.objective(ObjectiveSense.Minimize, Expr.dot(P,Q))

        if not quite:
            M.setLogHandler(sys.stdout)

        M.solve()
        status = M.getSolverDoubleInfo("intpntOptStatus")
        print("solver status " + str(status))
        if status == 1:
            lumped_time = M.getSolverDoubleInfo("intpntTime")
        else:
            lumped_time = 'inf'   


direct_scaling(np.eye(4),4,4)



def fixed_scaling(A,blk_size):
    pass
    
def search_for_scaling(A,P):
    pass

def ricaati_based(A,N,):
    pass
    # w/ fixed scaligns, solve for the continuous
    # scipy.linalg.solve_continuous_are(a, b, q, r)[source]
    # add error handling
    
    # A'X + XA - XBR^-1B'X+Q=0) directly using a Schur decomposition method.
    # t = Timer(lambda: linalg.solve_lyapunov(A,Q))
        # X = np.matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))


def read_data(file, key):
    pass

def write_data(key, value):
    pass





# testing_stuff()