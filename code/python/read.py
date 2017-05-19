from numpy import genfromtxt
import numpy as np
# import scipy as sp
from scipy import linalg

A = genfromtxt('counter.csv', delimiter=',')
X=    np.diag([0.4489, 0.3801,0.3798])
# print(P.size)
# P=P.reshape(3,3)
print(A)
# print("A eigen", np.linalg.eig(A)[0])
print(X)
# qwqew=linalg.solve_lyapunov(A, -np.eye(4))
# print('backed out',A*qwqew+qwqew*A.transpose())
# print('direct', np.linalg.eig(qwqew)[0])

# print("lyap eigen", np.linalg.eig(X)[0])
Q=np.dot(X,A)
print(Q)
print('paap eigen',np.linalg.eig(Q)[0])
A_syme=A
for i in range(3):
	diff=A_syme[i,i]
	for j in range(3):
		print (i,j,A[i,i]-A[i,j]/A[j,j]*A[j,i])
		# if (A[i,i]-A[i,j]/A[j,j]*A[j,i]>=0):
		# diff=diff-A_syme[i,j]/A_syme[j,j]*A_syme[j,i]
	# print (i, diff)
# print(A[0,0]-A[0,1]/(A[1,1])*A[1,0]-A[0,2]/(A[2,2])*A[2,0]-A[1,2]/(A[2,2])*A[2,1])
# print(A[0,0]-A[1,2]/(A[2,2])*A[2,1])s


# print(A[1,1]-A[1,0]/(A[0,0])*A[0,1]-A[1,2]/(A[2,2])*A[2,1])
# print(A[2,2]-A[2,0]/(A[0,0])*A[0,2]-A[2,1]/(A[1,1])*A[1,2])

# print(np.linalg.eig(my_data)[0])