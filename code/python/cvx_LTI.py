from cvxpy import *
import numpy as np
from numpy import genfromtxt
# A = np.random.randn(4,4)
# A=-np.eye(4)
A = genfromtxt('A.csv', delimiter=',')

X = Semidef(4)
# Y = Semidef(10)
constr1=[]
for i in range(4):
	for j in range(4):
		if i!=j:
			constr1.append(X[i,j]==0)
obj = Minimize(0)
constr1.append(X@A+A.transpose()@X<=-np.eye(4))
prob = Problem(obj, constr1)
prob.solve()
if X.value!=None:
	np.savetxt("A.csv", A, delimiter=",")
	np.savetxt("X.csv", X.value, delimiter=",")
	# for i in range(4):
	# 	for j in range(4):
	# 		print (A[i,i]-A[i,j]/A[j,j]*A[j,i])
	print(X.value)
	print()



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


