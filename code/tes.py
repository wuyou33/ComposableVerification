import sys
import mosek
from   mosek.fusion import *
import numpy as np

A=np.eye(10)
blk_size=5
i=2
P=[]
with Model("direct_scaling") as M: 
	# A=Matrix.sparse(10,10,) 
	for i in range(3):
		P.append(M.variable(Domain.inPSDCone(6)))
	E1=P[1].asExpr()
	# print()
	# print(Expr.neg(P1))

	# zero_block=Expr.reshape(Expr.zeros(36))
	zero_block=Expr.reshape(Expr.ones(0),6,0)
	# slices=Expr.ones(1)

	stacked=Expr.hstack(P[2],zero_block)
	print(stacked)
	# print(stacked.slice([0,0],[3,3]))
	# print(stacked)
	# print(stacked.slice(p[]))
	# # print(stacked.shape())
	# sliced=stacked.slice([1,1],[3,2])
	# print(sliced)
# 
	# M.constraint(zero_block.slice(4,5),Domain.equalsTo(slices))
# 
	# P.reshape(6,6)
	# A[0:5,0:5]=P
	# PP=Var.reshape(P1,36,1)

	# P=Var.vstack(PP,Var.hstack(zero_block,P2))

	# # A = Matrix.eye(6)
	# M.constraint(Expr.neg(Expr.add(Expr.mul(P,A),Expr.mul(A.transpose(),P))),Domain.inPSDCone())
	# M.solve()
# print(zero_block)

# print(A[(i-1)*blk_size:i*blk_size,:])