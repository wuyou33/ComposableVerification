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
# Filename:  riccati_solver.py
#
# Description:
#            This file contains the implementation of the discrete 
#            Riccati Equation Solver.
#
# References:
#            Direct solution algorithm taken from:
#              Laub, "A Schur Method for Solving Algebraic Riccati Equations."
#              U.S. Energy Research and Development Agency under contract 
#              ERDA-E(49-18)-2087.
#           
#           Iterative Techniques:
#           Simplistic Newton solver taken from:
#             Fabbender and Benner, "Initializing Newton's Method for Discrete-Time
#             Algebraic Riccati Equations Using the Butterfly SZ Algorithm." 
#             Proceedings of the 1999 IEEE International Symposium on Computer Aided 
#             Control System Design, Hawaii, USA, August 22-27, 1999.  pp. 70-74.
#           
#           Cyclic Reduction solver taken from:
#             Bini and Iannazzo, "A Cyclic Reduction Method for Solving Algebraic
#             Ricatti Equations." Technical Report, Dipartimento di Matematica, 
#             Universita di Pisa, 2005.
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
from    lyapunov_solver import LyapunovSolver

#
# Class that inherits the base Solver to implement specific Riccati solver
# attributes and methods associated to the SciPy and math libraries.
#
class RiccatiSolver:
    
    ### Solver's constants ###
    ITERATION_LIMIT     = 10000
    RICCATI_EPSILON     = 1.0E-5
    RELAXATION          = 1.0
    
    ### Solver's integration methods ###
    Iterative_Method, Direct_Method, Cyclic_Method = range(3)   
    
    lyap = None #LyapunovSolver()                 # Lyapunov solver. Used in the Newton's iterative step.
    
    A = None #numpy.empty()                       # A is an input numpy Matrix
    B = None #numpy.empty()                       # B is an input numpy Matrix
    Q = None #numpy.empty()                       # Q is an input numpy Matrix
    R = None #numpy.empty()                       # R is an input numpy Matrix
    Integration_Method = Iterative_Method   # Default Riccati method used by the solver
    iteration_limit    = ITERATION_LIMIT    # Default iteration limit when using iterative method
    use_cyclic         = False              # Use cyclic reduction in interactive method?
    x                  = 0                  #
    dx                 = 0                  # 
    a0   = None #numpy.empty()  
    a1   = None #numpy.empty()  
    h    = None #numpy.empty()
    k    = None #numpy.empty()
    hhat = None #numpy.empty()     
      
       
    #
    # Init(A, B, Q, R, method = Iterative_Method, iteration_limit = ITERATION_LIMIT)
    # Initialization/setup method
    # $Inputs:
    #    A:                Input matrix. 
    #    B:                Input matrix.
    #    Q:                Input Matrix
    #    R:                Input Matrix    
    #    method:           Selected method used to solve the Riccati equations
    #    use_cyclic:       Use cyclic reduction instead Newton's method during the interactive method?
    #                      In case of a Direct method, it is used the Newton method.
    #
    # $Outputs:
    #    none
    #
    def Init(self, A, B, Q, R, method = Iterative_Method, use_cyclic = False, iteration_limit = ITERATION_LIMIT):       
        self.A = A
        self.B = B
        self.Q = Q
        self.R = R
        self.Integration_Method = method
        self.iteration_limit = iteration_limit
        self.use_cyclic = use_cyclic        
        self.iterations = 0  
        self.lyap = LyapunovSolver(A, Q) 
    
    #
    def __init__(self, A, B, Q, R, method = Iterative_Method, use_cyclic = False, iteration_limit = ITERATION_LIMIT):
        self.Init(A, B, Q, R, method, use_cyclic, iteration_limit) 
    
    #
    # DirectSolver()
    # It Solves the Discrete Riccati Equation (A'XA-(A'XB)(R+B'XB)^-1(B'XA) -X +Q = 0) using a direct
    # approach.  The function returns an estimative of X based on the A and Q input
    # matrices.  
    # It is used the Schur decomposition method.
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    #        
    def DirectSolver(self):
        
        g   = self.B * numpy.linalg.inv(self.R) * self.B.transpose()
        fit = numpy.linalg.inv(self.A).transpose()
        
        z11 = self.A + g * fit * self.Q
        z12 = -1.0 * g * fit
        z21 = -1.0 * fit * self.Q
        z22 = fit
        z   = numpy.vstack((numpy.hstack((z11, z12)), numpy.hstack((z21, z22))))
        
        # s:    Schur form of A. It is real-valued for the real Schur decomposition.
        # u:    An unitary Schur transformation matrix for A.
        # sdim: Dimension of s.
        # lhp:  Left Hand Plane (x < 0)        
        [s,u,sdim] = scipy.linalg.schur(numpy.linalg.inv(z), sort='lhp')
        
        (m,n) = u.shape
        
        u11 = u[0:m/2, 0:n/2]
        u12 = u[0:m/2, n/2:n]
        u21 = u[m/2:m, 0:n/2]
        u22 = u[m/2:m, n/2:n]
        u11i = numpy.linalg.inv(u11)
        
        solution =  numpy.asmatrix(u21)*numpy.asmatrix(u11i)
        
        return solution
    
    
    #
    # IterativeSolverCyclicStep()
    # Steps the cyclic reduction solver one iteration. 
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    #
    def IterativeSolverCyclicStep(self):
        hinv = numpy.linalg.inv(self.h)
        h1 = self.h - \
             self.k*hinv*self.k.transpose() - \
             self.k.transpose()*hinv*self.k
             
        self.hhat = self.hhat - self.k*hinv*self.k.transpose()
        
        self.k = -1.0*self.k*hinv*self.k
        
        self.h = h1 
    
    
    #
    # IterativeSolverCyclic()
    # Initializes the cyclic reduction solver variables. 
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    #
    def IterativeSolverCyclicInit(self):
        
        binv    = numpy.linalg.inv(self.B)
        
        self.a1 = -1.0 * self.A.transpose() * binv.transpose() * self.R * binv
        
        self.a0 = binv.transpose() * self.R * binv + \
                  self.A.transpose() * binv.transpose() * self.R * binv * self.A + \
                  self.Q
                  
        self.h = self.a0
        self.k = self.a1
        self.hhat = self.h
        
    #
    # IterativeSolverCyclic()
    # Solves the Riccati equation using the cyclic reduction method.
    # B Matrix must be square. Depending on the problem,
    # this method could be faster than Newton interactive method. 
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    #     
    def IterativeSolverCyclic(self):
               
        self.IterativeSolverCyclicInit()
        
        error = 1.0E+6
        count = 0
        while (error > self.eps and count < self.iteration_limit) or count < 2:
            
            self.IterativeSolverCyclicStep()
            error = numpy.linalg.norm(self.k)
            count = count + 1
        
        z = -1.0 * numpy.linalg.inv(self.hhat) * self.a1
        binv = numpy.linalg.inv(self.B)
        try:
            zinv = numpy.linalg.inv(z)
        except:
            warnings.warn('During the cyclic reduction, a singular matrix was found - using psuedo-inverse', RuntimeWarning)
            zinv = numpy.linalg.pinv(z) 
        
        solution = binv.transpose()*(self.R * binv * (self.A - z)) * zinv
        
        return solution
    
    #
    # IterativeSolverNewtonCost(x)
    # Computes the current error in the Riccati solution estimate for use with
    # the Newton iterative solver
    # $Inputs:
    #    x:    number of steps
    # $Outputs:
    #    Matrix of costs
    #
    def IterativeSolverNewtonCost(self, x):       
        solution = self.Q - x + self.A.transpose() * x * self.A - \
                   self.A.transpose() * x * self.B * numpy.linalg.inv(self.R + self.B.transpose() \
                   * x * self.B) * self.B.transpose() * x * self.A
        
        return solution

    #
    # IterativeSolverNewtonStep()
    # Steps the Newton iterative solver one iteration.
    # $Inputs:
    #    none
    # $Outputs:
    #    The number of steps/iterations
    #    
    def IterativeSolverNewtonStep(self):
        ak = self.A - self.B * numpy.linalg.inv(self.R + self.B.transpose() * self.x * self.B) \
             * self.B.transpose() * self.x * self.A
        
        # The iterative Lyapunov solver must be used here due to the possible presence of
        # numerical instabilities.  However, the necessary accuracy at this step can be
        # low, so the iteration count is set to 30.
        self.lyap.Init(ak.transpose(), self.IterativeSolverNewtonCost(self.x))  
        self.dx = self.lyap.Solve()
        
        self.x = self.x + self.RELAXATION * self.dx
    
    #
    # IterativeSolverNewton()
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    # 
    def IterativeSolverNewton(self):
        self.x = numpy.eye(self.Q.shape[0])
        error = 1.0E+6
        iterations = 0
        while error > self.RICCATI_EPSILON and iterations < self.iteration_limit:
            
            self.IterativeSolverNewtonStep()
            error = abs(self.dx.max())
            iterations = iterations + 1
        
        solution = self.x
        
        return solution 
    
          
    #
    # IterativeSolver()
    # It Solves the Discrete Riccati Equation (A'XA-(A'XB)(R+B'XB)^-1(B'XA) -X +Q = 0) using an iterative
    # approach.  The function returns an estimative of X based on the A and Q input
    # matrices.  
    # This iterative solver requires that the eigenvalues of the square matrix
    # A be within the unit circle for convergence reasons.
    # $Inputs:
    #    none
    # $Outputs:
    #    none
    #        
    def IterativeSolver(self, use_cyclic = False): 
        
        if self.use_cyclic == False:
            return self.IterativeSolverNewton()
        else:
            (m,n) = self.B.shape
            if m != n:
                warnings.warn("For Cyclic reduction method, 'B' matrix must be square - Using Newton method", RuntimeWarning)
                return self.IterativeSolverNewton()
            else:
                return self.IterativeSolverCyclic()
    
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
        