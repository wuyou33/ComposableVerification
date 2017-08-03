# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Filename:  nonlinear_model.py
#
# Description:
#            This file sets up a symbolic math based nonlinear model
#            as well as tools to convert it to a linearized form.
#
# Updates:
#            Date/Time        Author                 Description
#            Oct 26, 2016     Gray Thomas            File was created
#

import cvxpy as cvx
import numpy as np
import control as ctrl
import sympy as sp

# constants from the paper TedrakeManchesterTobenkinRoberts2010
g=9.8 # meters/second/second
l=0.5 # meters
b=0.1 # meters*meters*Kg/second
m=1.0 # Kg
Q = np.diagflat([10.0,1.0])
R = np.eye(1)*15.0

print "First we define some symbolic varables"
x = sp.Symbol("x")
v = sp.Symbol("v")
u = sp.Symbol("u")
X = np.array([[x], [v]])
# pendulum dynamics: deriv[v] = -gl sin(x) 
# An affine nonlinear system:
# ddot{x} = f(x)+g(x)u 
print "we use them to define the nonlinear affine form:"
f = np.array([[v],[-b*v-g*l*sp.sin(x)]])
print "f(x,v) = :\n",f
g = np.array([[sp.Float(0)],[sp.Float(1/m/l/l)]])
print "g(x,v) = :\n",g
print 

def factorial(n):
    if n<=1:
        return 1
    return n*factorial(n-1)

def part_taylor(func, X, nums, X0):
    return sp.N(func.diff(X[0,0],nums[0]).subs(X[0][0],X0[0,0]).diff(X[1,0],nums[1]).subs(X[1][0],X0[1,0])/factorial(nums[0])/factorial(nums[1]))

def sparse_taylor(func, var, var0, N=10, tol=1e-9):
    ret=[]
    for S in range(N+1):
        for i in range(S+1):
            j=S-i
            coeff = part_taylor(func,var,(i,j),var0)
            if abs(coeff)>tol:
                ret.append((i,j,coeff))
    return ret

def subs(obj, X, X0):
    for i in range(X.shape[0]):
        obj = obj.subs(X[i,0],X0[i,0])
    return float(sp.N(obj))

def lin(f, X, X0):
    A = np.array([
        [float(part_taylor(f[0,0], X, (1,0), X0)),float(part_taylor(f[0,0], X, (0,1), X0))],
        [float(part_taylor(f[1,0], X, (1,0), X0)),float(part_taylor(f[1,0], X, (0,1), X0))]])
    B = np.array([
        [subs(g[0,0], X, X0)],
        [subs(g[1,0], X, X0)]])
    return A, B