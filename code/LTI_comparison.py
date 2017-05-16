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

        # Setting up the constraints, i is the constraint index
        for i in range(blks):
            row_slice = Expr.neg(Expr.mul(P[i],A[(i)*blk_size:(i+1)*blk_size,:]))
            colum_slice = row_slice.transpose()
            if i==0:
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
            # scalings = construct_blk_diag(N,blk_size,i)
            PAAP=(Expr.add(rowstaced,columnstaced))
        # adding the scaling portion

        # for i in range(blks):
        #     zero_on_left=Expr.reshape(Expr.zeros(),,)
        #     zero_on_right=Expr.reshape(Expr.zeros(),,)    
        #     for j in range(blks):            
        #         if j==0 and i==0:
        #             # adding an neg scaling to 
        #             neg_scaling=Expr.neg(N[i])
        #             neg_scaling_row=Expr.hstack(zero_on_left,neg_scaling,zero_on_right)
        #             row_slice=Expr.add(row_slice,neg_scaling)
        #         else:
        #             # adding the original positive scaling to 
        #             pos_scaling_row=Expr.hstack(zero_on_left,positive,zero_on_right)
        #             row_slice=Expr.add(row_slice,pos_scaling_row)
        #             next_one = Expr.hstack(this_one, next_one)
        #             continue


        # setting up the constraints, essentially solving for n independent LMIs
            M.constraint(PAAP,Domain.inPSDCone())

        if not feasibility_only:
            M.objective(ObjectiveSense.Minimize, Expr.dot(P,Q))

        if not quite:
            M.setLogHandler(sys.stdout)

        M.solve()
        status = M.getSolverDoubleInfo("intpntOptStatus")
        print("solver status " + str(status))
        if status == 1:
            lumped_time = M.getSolverDoubleInfo("intpntTime")

            # print("solver lumped_time " + str(lumped_time))
            # print(P[2].level())
        else:
            lumped_time = 'inf'   


# direct_scaling(np.random.randn(6,6),6,2)

def construct_blk_diag(N, blks, blk_size, index):
    new_blks=[]
    zero_blks = Expr.reshape(Expr.zeros((blks)*blk_size*blk_size),blk_size,(blks)*blk_size)
    for i in range(blks):
        # print(i)
        # print(zero_blks)
        # print(N[i+1])
        if i ==0:
            new_blks.append(N[i])
        else:
            new_blks.append(Expr.hstack(zero_blks,N[i]))

    final=N[0].asExpr()
    for i in range(blks-1):
        # print(i)
        final=Expr.hstack(final,new_blks[i+1])
        # print(final)
    # final=Expr.hstack(final,new_blks[blks-1])

    print('before reshaping')
    # print(final.slice([0,0],[2,26]))

    # print(oneseew.slice([0,0],[1,1]))

    # final=Expr.reshape(final,blks*blk_size,blks*blk_size)
    # print("after reshaping")
    # print(final)
    return final


    # for i in range(blks):
    #     if index==0:
    #         first = Expr.neg(N[0])
    #     else:
    #         first = N[0]


    #     if i == index:
    #         incremental_block=Expr.hstack(Expr.neg(N[i]),zero_blks)
    #     else:
    #         incremental_block=Expr.hstack(N[i],zero_blks)
    #     blk_diag_form=Expr.hstack(blk_diag_form,incremental_block)
    #     i+=1

    blk_diag_form=blk_diag_form.reshape(blks*blk_size,blks*blk_size)
    return blk_diag_form 

blks=3
blk_size=3
N=[]
with Model("testing") as M:    
    for i in range(blks):
        N.append(M.variable(Domain.inPSDCone(blk_size)))
    final = construct_blk_diag(N,blks,blk_size,0)
    ones=Matrix.eye(blks*blk_size)



def fixed_scaling(A,blk_size):
    pass
    
def search_for_scaling(A,P,N,):
    pass

def ricaati_based(A,N,):
    pass
    # w/ fixed scaligns, solve for the continuous
    # scipy.linalg.solve_continuous_are(a, b, q, r)[source]
    # add error handling
    
    # t = Timer(lambda: linalg.solve_lyapunov(A,Q))
        # X = np.matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))


def read_data(file, key):
    pass

def write_data(key, value):
    pass





# testing_stuff()