import numpy as np
import scipy.linalg


def comp():
# A’X + XA - XBR^-1B’X+Q=0

    #first, try to solve the ricatti equation
    n = np.size(A)
    Q = np.identity(n)
    R = np.identity(n)
    X = np.matrix(scipy.linalg.solve_continuous_are(A, B, Q, R))
     
    #compute the LQR gain
    K = np.matrix(scipy.linalg.inv(R)*(B.T*X))
     
    eigVals, eigVecs = scipy.linalg.eig(A-B*K)
     
    return K, X, eigVals
