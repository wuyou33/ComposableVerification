# Decision Making Framework Simulation_Based on LQR-Tree
# Copyright (C) 2016
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Filename:  shooting_sequential_convex_optimization.py
#
# Description:
#            This file contains a starting implementation of a discrete
#            time input vector shooting with sequential convex MPC.
#
# Updates:
#            Date/Time        Author                 Description
#            Oct 26, 2016     Gray Thomas            File was created
#

# Importing specific modules
import cvxpy as cvx
import numpy as np
import time
import control as ctrl
import sympy as sp
from sympy.mpmath import taylor
from math import pi, exp
from nonlinear_model import f, g, X, u, subs, lin
import matplotlib.pyplot as plt

def integrate_model(A,B,dt,subn=100):
    """
    integrate_model(A,B,dt,subn=100):
    xdot = A x + B u
    becomes 
    x_{k+1}= Aint x_k + Bint u_k
    via approximate matrix exponentiation

    return Aint, Bint
    """
    xdim=A.shape[0]
    udim=B.shape[1]
    Aint=np.eye(xdim)
    Bint=np.zeros((xdim,udim))
    for i in range(0,subn):
        Aint+=(1.0/subn)*dt*A.dot(Aint)
        Bint+=(1.0/subn)*dt*Aint.dot(B)
    return Aint, Bint

tmax = 2.0
steps = 10
dt = tmax / (steps)
assert abs(sum([dt] * steps) - tmax) < 1e-7


def check_inclusion(state,rho=10.0):
    S = np.array([[ 26.73559192 , 10.11439662],
        [ 10.11439662,   4.36840001]])
    return (state-np.array([[np.pi],[0]])).T.dot(S).dot(state-np.array([[np.pi],[0]]))<rho


def sim(state0, input_vector):
    lins = []
    Xs = [state0]
    Us = []
    state = np.array(state0)
    for i in range(steps):
        A, B = lin(f, X, state)
        xdot = subs(f[0, 0], X, state) + subs(g[0, 0],
                                              X, state) * input_vector[i]
        vdot = subs(f[1, 0], X, state) + subs(g[1, 0],
                                              X, state) * input_vector[i]
        state0 = np.array(state)
        state = state + dt * np.array([[xdot, vdot]]).T
        xdotF = subs(f[0, 0], X, state) + subs(g[0, 0],
                                               X, state) * input_vector[i]
        vdotF = subs(f[1, 0], X, state) + subs(g[1, 0],
                                               X, state) * input_vector[i]
        AF, BF = lin(f, X, state)
        xdotavg = 0.5 * (xdot + xdotF)
        vdotavg = 0.5 * (vdot + vdotF)
        Aavg = 0.5 * (A + AF)
        Bavg = 0.5 * (B + BF)
        state = state0 + dt * np.array([[xdotavg, vdotavg]]).T
        Adisc, Bdisc = integrate_model(Aavg, Bavg, dt, subn=5)
        lins.append((Adisc, Bdisc))
        Xs.append(state)
        Us.append(input_vector[i])
        if check_inclusion(state, rho=1) and i>1:
            print "truncated at ", len(Us)
            print state.T
            break
    return Xs, Us, lins



def cost(Xs, input_vector, goal_state):
    cost = 0.0
    at_goal = True
    for i in reversed(range(steps)):
        if not np.allclose(Xs[i], goal_state, rtol=1e-4):
            at_goal = False
        if not at_goal:
            cost += input_vector[i]**2 * dt
            cost += 1.0 * dt
    return cost




def convex_opt(Xs, Us, lins, max_x_dev, max_u):
    dX = cvx.Variable(2, len(Us) + 1)
    dU = cvx.Variable(1, len(Us))
    constraints = []
    Q=cvx.Constant(np.array([[1.0]]))
    R = cvx.Constant(np.eye(2))
    obj = 0.0
    for i in range(len(Us)):
        A, B = lins[i]

        constraints.append(dX[:, i + 1]
                           == cvx.Constant(A) * dX[:, i] 
                           + cvx.Constant(B) * dU[:, i])
        constraints.append(Us[i] + dU[0, i] > -max_u)
        constraints.append(Us[i] + dU[0, i] < max_u)
        constraints.append(dX[0, i] > - max_x_dev)
        constraints.append(dX[0, i] < max_x_dev)
        obj += cvx.abs(Us[i]+dU[0,i])
        # obj += cvx.abs(Xs[i][0,0]+dX[0,i+1] - goal_X[0,0])*cvx.Constant(lam*exp(i*dt*omega))
        # obj += cvx.abs(Xs[i][1,0]+dX[1,i+1] - goal_X[1,0])*cvx.Constant(lam*exp(i*dt*omega))
    constraints.append(dX[0, 0] == 0.0)
    constraints.append(dX[1, 0] == 0.0)
    # constraints.append(dX[0, len(Us)] == 0.0)
    # constraints.append(dX[1, len(Us)] == 0.0)
    goal=np.array([[np.pi],[0.0]])
    S = np.array([[ 26.73559192 , 10.11439662],
        [ 10.11439662,   4.36840001]])
    constraints.append(cvx.quad_form(Xs[len(Us)]+dX[:,len(Us)]-goal,S) <= cvx.Constant((Xs[len(Us)]-goal).T.dot(S).dot((Xs[len(Us)]-goal))))
    obj += cvx.quad_form(Xs[len(Us)]+dX[:,len(Us)]-goal,S)*10.0
    prob=cvx.Problem(cvx.Minimize(obj), constraints)
    prob.solve(solver=cvx.MOSEK, verbose=True)
    X_star = [Xs[i]+dX.value[:,i].reshape((2,1)) for i in range(len(Us)+1)]
    U_star = [Us[i]+dU.value[0,i] for i in range(len(Us))]
    print dX.value
    # plt.plot(dX.value.T)
    # plt.plot([X_star[i][0,0] for i in range(steps+1)])
    # plt.show()
    
    return X_star, U_star


goal_state = np.array([[np.pi, 0.0]]).T
# discrete_input_vector = [2.0] * 10 + [-2.0]*10+[0.0]*( steps - 20)


max_x_dev = 0.5  # index in x, constraint on dX
max_u = 2.0
for j in range(3):
    for tries in range(800):
        state0 = np.array([[np.random.normal(1.0), np.random.normal(2.0)]]).T
        while check_inclusion(state0, rho=40) or not check_inclusion(state0, rho=90):
            state0 = np.array([[np.random.normal(1.0), np.random.normal(2.0)]]).T
        assert not check_inclusion(state0)
        Us0 = [np.random.uniform(-2,2) for i in range(steps)]
        Xs_star, Us_star, LL = sim(state0, Us0)
        break

    # if len(Us_star)==steps:
    #     raise Exception()
    fig, axs = plt.subplots(2,1,sharex=True)
    Xtrajectory_history=[Xs_star]
    Utrajectory_history=[Us_star]
    linearization_history = [LL]
    for i in range(10):
        # plt.plot([x for x in Us_star])
        Xs_starp, Us_starp = convex_opt(Xtrajectory_history[-1], Utrajectory_history[-1],linearization_history[-1], max_x_dev, max_u)
        Xs_starpp, Us_starpp, LLp = sim(state0, Us_starp)
        Utrajectory_history.append(Us_starpp)
        Xtrajectory_history.append(Xs_starpp)
        linearization_history.append(LLp)

    for Xs_star, Us_star, lins in zip(Xtrajectory_history,Utrajectory_history,linearization_history):
        axs[0].plot([x for x in Us_star])
        axs[1].plot([x[0,0] for x in Xs_star])



# print cost(Xs, discrete_input_vector, goal_state)
plt.show()
