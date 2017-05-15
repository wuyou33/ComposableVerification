import cvxpy as cvx
import numpy as np
W = cvx.Variable(3,3)
S= cvx.Variable(3,3)
B = np.zeros([3,6])
# T = cvx.diag(cvx.vstack(B, W))
T=cvx.hstack(W,S,B)
print(T)