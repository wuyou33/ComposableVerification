# Decision Making Framework based on Social Cognition Theory (SCT)
# Copyright (C) 2016 Hilgad Montelo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Filename:  lyapunov_solver.py
#
# Description:
#            This file contains the implementation of the discrete 
#            Lyapunov Equation Solver.
#
# References:
#           Davinson and Man, "The Numerical Solution of A'Q+QA=-C." 
#            IEEE Transactions on Automatic Control, Volume 13, Issue 4, August, 1968.  p. 448.
#
#           Discrete-Time Lyapunov Solver based on Octave dlyap.m file
#             Copyright (C) 1993, 1994, 1995, 2000, 2002, 2004, 2005, 2007
#             Auburn University.  All rights reserved.
#           
#           Uses Schur decomposition method as in Kitagawa,
#             "An Algorithm for Solving the Matrix Equation @math{X = F X F' + S}}",
#             International Journal of Control, Volume 25, Number 5, pages 745--753
#             (1977).
#           
#           Column-by-column solution method as suggested in
#             Hammarling, "Numerical Solution of the Stable, Non-Negative
#             Definite Lyapunov Equation", Journal of Numerical Analysis, Volume
#             2, pages 303--323 (1982).
#
#            Jeffrey Armstrong: Approximatrix, LLC
#
# Updates:
#            Date/Time        Author                 Description
#            Oct 21, 2016     Hilgad.Montelo         File was created
#

# Importing specific modules
import  numpy.linalg
import  scipy.linalg
import  warnings

#
# Class that inherits the base Solver to implement specific Lyapunov solver
# attributes and methods associated to the SciPy and math libraries.
#
class LyapunovSolver:
    
    ### Solver's constants ###
    ITERATION_LIMIT     = 10000
    LYAPUNOV_EPSILON    = 1.0E-6
    
    ### Solver's integration methods ###
    Iterative_Method, Direct_Method = range(2)
    
    A = None #numpy.empty()                       # A is an input Matrix
    Q = None #numpy.empty()                       # Q is an input Matrix
    Integration_Method = Iterative_Method   # Default Lyapunov method used by the solver
    iteration_limit    = ITERATION_LIMIT    # Default iteration limit when using iterative method
    eps                = LYAPUNOV_EPSILON
    
    
    #
    # Init(a,q,method=Direct_Method,iteration_limit=ITERATION_LIMIT)
    # Initialization/setup method
    # $Inputs:
    #    A:                Input matrix. 
    #    Q:                Input Matrix
    #    method:           Selected method used to solve the Lyapunov equations
    #    iteration_limit:  Iteration limit used in the Lyapunov's iterative method
    #
    # $Outputs:
    #    none
    #
    def Init(self, A, Q, method = Iterative_Method, iteration_limit = ITERATION_LIMIT, eps = LYAPUNOV_EPSILON):
        self.A = A
        self.Q = Q
        self.Integration_Method = method
        self.iteration_limit = iteration_limit
        self.eps = eps
        
    #
    def __init__(self, A, Q, method = Iterative_Method, iteration_limit = ITERATION_LIMIT, eps = LYAPUNOV_EPSILON):
        self.Init(A, Q, method, iteration_limit, eps)
    
    #
    # DirectSolver()
    # It Solves the Discrete Lyapunov Equation (A X A' -X + Q = 0) using a direct
    # approach.  The function returns an estimative of X based on the A and Q input
    # matrices.  
    # It is used the Schur decomposition method.
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    #        
    def DirectSolver(self):
        (m,n) = self.A.shape
        if m != n:
            raise ValueError("Matrix 'A' must be a square matrix") 
        
        # s:    Schur form of A. It is real-valued for the real Schur decomposition.
        # u:    An unitary Schur transformation matrix for A.
        # sdim: Dimension of s.
        # lhp:  Left Hand Plane (x < 0)
        [s, u, sdim] = scipy.linalg.schur(self.A, sort='lhp')
        
        s = numpy.asmatrix(s)
        u = numpy.asmatrix(u)
        
        b = u.transpose()*self.Q*u
        
        x = numpy.asmatrix(numpy.zeros(self.A.shape))
        
        j = n-1
        
        while j >= 0:
            
            j1 = j+1
            
            # Verifies the Schur block
            blocksize = 1
            if j > 0:
                if s[j,j-1] != 0.0:
                    blocksize = 2
                    j = j - 1
            
            Ajj = scipy.linalg.kron(s[j:j1,:][:,j:j1] ,s) - numpy.eye(blocksize*n)
            
            rhs = numpy.reshape(b[:,j:j1], (blocksize*n, 1))
            
            if j1 < n:
                rhs2 = s*(x[:,j1:]*(s[j:j1,:][:,j1:]).transpose())
                rhs = rhs + numpy.reshape(rhs2,(blocksize*n, 1))
                
            v = -1.0 * scipy.linalg.solve(Ajj,rhs)
            
            x[:,j] = v[0:n]
            
            if blocksize == 2:
                x[:,j1-1] = v[n:blocksize*n+1]
            
            j = j - 1
        
        x = u * x * u.transpose()
        
        return x
       
    #
    # IterativeSolver()
    # It Solves the Discrete Lyapunov Equation (A X A' -X + Q = 0) using an iterative
    # approach.  The function returns an estimative of X based on the A and Q input
    # matrices.  
    # This iterative solver requires that the eigenvalues of the square matrix
    # A be within the unit circle for convergence reasons.
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    #        
    def IterativeSolver(self):
        
        error = 1E+6
        
        x = self.Q
        ap = self.A
        apt = self.A.transpose()
        count = 1
        
        (m,n) = self.A.shape
        if m != n:
            raise ValueError("Matrix 'A' must be a square matrix") 
        
        det_a = numpy.linalg.det(self.A)
        if det_a > 1.0:
            raise ValueError("Matrix 'A' must have eigenvalues within the unit circle") 
        
        while error > self.eps and count < self.iteration_limit:
            
            change = ap*self.Q*apt                
            x = x + change
            
            ap = ap*self.A
            apt = apt*(self.A.transpose())
            error = abs(change.max())
            count = count + 1
        
        if count >= self.iteration_limit:
            warnings.warn("Lyapunov Iterative Solver: Iterative Limit Reached - No convergence", RuntimeWarning)
            
        return x  
    
    #
    # Solve()
    # Method used to calculate the isystem's equations.
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    #
    def Solve(self):        
        if self.Integration_Method == self.Iterative_Method:
            return self.IterativeSolver()      
        else:
            return self.DirectSolver()
        