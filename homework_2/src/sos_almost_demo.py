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
# Filename:  sos_almost_demo.py
#
# Description:
#            This file contains a starting implementation of the algorithm presented
#            in the following paper:
#            Russ Tedrake, Ian R. Manchester, Mark Tobenkin, John W. Roberts, "LQR-Trees: Feedback Motion
#              Planning viaSums-of-Squares Verification." Computer Science and Artificial Intelligence Lab.
#              Massachusetts Institute of Technology. February 27, 2010.
#
# Updates:
#            Date/Time        Author                 Description
#            Oct 18, 2016     Gray Thomas            File was created
#            Oct 22, 2016     H. Montelo             Updated File
#

# Importing specific modules
import cvxpy as cvx
import numpy as np
import time
import control as ctrl
import sympy as sp
from sympy.mpmath import taylor
from math import pi

from nonlinear_model import *

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from math import sin, cos, sqrt
from riccati_solver import RiccatiSolver


# def part_taylor(func, x, numx, y, numy):
#     return func.diff(x,numx).diff(y,numy)/factorial(numx+numy)


def string_from_dict(dd, names=["x", "y", "z"], precision=2):
    ret = []
    fmt = "(%%.%de)" % precision
    keys = dd.keys()
    for i in range(len(list(dd.keys()[0]))):
        keys = sorted(keys, key=lambda tup: tup[i])
    keys = sorted(keys, key=lambda tup: sum(tup))
    for monom in keys:
        coeff = dd[monom]
        stringit = ""
        if isinstance(coeff, float):
            if abs(coeff) < 1e-9:
                continue
            stringit += fmt % coeff
        else:
            # try:
            stringit += "'" + coeff.name() + "'"
            # except :
            #     pass
            # stringit += repr(type(coeff))+repr(coeff)
        for j, q in enumerate(list(monom)):
            if q == 1:
                stringit += names[j]
            elif q > 1:
                stringit += names[j] + "^%d" % q
        ret.append(stringit)
    return " + ".join(ret)

# print dir(cvx.Variable()+cvx.Variable("h2"))
# exit()


def string_from_taylor(tay, precision=2):
    ret = []
    fmt = "(%%.%de)" % precision
    for i, j, coeff in tay:
        sx, sy, scof = "", "", ""
        if abs(coeff) > 1e-9:
            scof = fmt % coeff
            if i == 1:
                sx = "dx"
            elif i > 1:
                sx = "dx^%d" % i
            if j == 1:
                sy = "dv"
            elif j > 1:
                sy = "dv^%d" % j
            ret.append(scof + sx + sy)
    return " + ".join(ret)


def PrintUTHeader():
    print "University of Texas at Austin"
    print "Mechanical Engineering Department"
    print "ME_396D - Decision_Control-Human_Centered_Robotics"
    print "Fall 2016, Homework #2"
    print "Hilgad Montelo"
    print "Gray Thomas"
    print " "
    print "--- LQR-Trees: Feedback Motion Planning via SOS Verification ---"
    print " "


def add_tuple(t1, t2):
    return tuple(x1 + x2 for x1, x2 in zip(t1, t2))


def poly_mult(poly1, poly2):
    k1 = poly1.keys()
    k2 = poly2.keys()
    poly3 = dict()
    for m1 in k1:
        for m2 in k2:
            m3 = add_tuple(m1, m2)
            if not poly3.has_key(m3):
                poly3[m3] = poly1[m1] * poly2[m2]
            else:
                poly3[m3] += poly1[m1] * poly2[m2]
    return poly3


def poly_add(poly1, poly2):
    new_poly = dict(poly1)
    for k in poly2.keys():
        if new_poly.has_key(k):
            new_poly[k] += poly2[k]
        else:
            new_poly[k] = poly2[k]
    return new_poly


def listtodict(lst):
    dct = {}
    for q in lst:
        dct[q[:-1]] = q[-1]
    return dct


def nosympy(dct):
    for key in dct.keys():
        dct[key] = float(dct[key])
    return dct


def get_monomial_lookup(monomials):
    dct = {}
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            monom = add_tuple(monomials[i], monomials[j])
            if not dct.has_key(monom):
                dct[monom] = []
            dct[monom].append((i, j))
    return dct


def setupSOSconstraint(Matrix, poly_dct, monomials):
    lookup = get_monomial_lookup(monomials)
    constraints = []
    for monom in lookup.keys():
        if poly_dct.has_key(monom):
            coeff = poly_dct[monom]
        else:
            coeff = 0.0
        constraints.append(coeff == sum(
            [Matrix[i, j] for i, j in lookup[monom]]))
    return constraints


def plotrho(rhos, color='k', order=4):
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, LightSource
    from matplotlib.ticker import MaxNLocator
    import matplotlib.cm as cm
    from math import sin, cos, sqrt
    import numpy as np

    # make these smaller to increase the resolution
    Qdx, Qdy = 0.2, 0.4

    x_range = [-2, 6]
    y_range = [-13, 13]

    # generate 2 2d grids for the x & y bounds
    Qy, Qx = np.mgrid[slice(y_range[0], y_range[1] + Qdy, Qdy),
                      slice(x_range[0], x_range[1] + Qdx, Qdx)]
    Qz = Qx * 0
    Qz2 = Qx * 0

    dVp = Dpoly(nosympy(listtodict(sparse_taylor(dV, dX, dX0, N=order, tol=1e-9))))
    print dVp
    dVpk = dVp.keys()
    # calculate analytic second order approximation to Vdot


    print A
    print B
    print Q
    print R
    print P
    Z=P.dot(B).dot(np.linalg.inv(R)).dot(B.T).dot(P) + Q
    print Z*10.11439662/47.12931231

    e, U = np.linalg.eigh(P)
    print -Z
    print U.dot(np.diagflat(e)).dot(U.T)
    print U
    print np.array([[-.9,-.3]]).dot(Z).dot(np.array([[-.9],[-.3]]))
    print np.array([[-.3,.9]]).dot(Z).dot(np.array([[-.3],[.9]]))
    from math import atan2
    print atan2(.3,.9)*180/np.pi
    print 
    print e
    print U.T
    print U.dot(U.T)

    print dVp
    print Z
    print "Z matches dVp"



    for r in range(Qx.shape[0]):
        for c in range(Qx.shape[1]):
            numeric_dx = Qx[r, c] - X0[0, 0]
            numeric_dv = Qy[r, c] - X0[1, 0]
            dVp(numeric_dx, numeric_dv)
            Qz[r,c]=float(dV.subs(dx,numeric_dx).subs(dv, numeric_dv))
            # Qz[r, c] = float(dV.subs(dx, numeric_dx).subs(dv, numeric_dv))
            # Qz2[r, c] = float(dV.subs(dx, numeric_dx).subs(dv, numeric_dv))
            Qz2[r, c] = dVp(numeric_dx, numeric_dv)

            # Qz[r,c] = -np.array([[numeric_dx, numeric_dv]]).dot(Z).dot(np.array([[numeric_dx], [numeric_dv]]))
            # for second order should equal dx^T P dx
            # Qz[r,c]=Qz[r,c]
            # print Qz[r,c]

    Pprime = P

    KP = np.linalg.cholesky(Pprime)

    assert np.linalg.norm(Pprime-KP.dot(KP.T))<1e-7
    

    # x and y are bounds, so z should be the value *inside* those bounds.
    # Therefore, remove the last value from the z array.
    Qz = Qz[:-1, :-1]
    Qz2 = Qz2[:-1, :-1]

    levels = MaxNLocator(nbins=15).tick_values(Qz.min(), Qz.max())

    # pick the desired colormap, sensible levels, and define a normalization
    # instance which takes data values and translates those into levels.
    cmap = plt.get_cmap('PiYG')
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    cmap.set_over("r", 1.0)
    cmap.set_under("g", 1.0)

    fig, (ax0, ax1) = plt.subplots(nrows=2)
    ve = 1
    ls = LightSource(azdeg=315, altdeg=25)
    rgb = ls.shade(Qz, cmap=cmap)
    # im = ax0.pcolormesh(Qx, Qy, Qz, cmap=cmap, norm=norm)
    im = ax0.imshow(Qz, interpolation='bilinear', cmap=cmap,
                    origin='lower', extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
                    vmax=0.0, vmin=Qz.min(), aspect='auto')
    # im = ax0.contour(Qx[:-1, :-1] + Qdx/2.,
    #                   Qy[:-1, :-1] + Qdy/2., Qz,
    #                   cmap=cmap)
    fig.colorbar(im, ax=ax0)
    ax0.set_title(r'$\hat{\dot{V}}(x,v)$')
    for rho in rhos:
        iso_x, iso_y = [], []
        for o in np.linspace(0, 2 * np.pi, 1000):
            p = np.linalg.solve(KP.T, np.array([[cos(o), sin(o)]]).T)
            p *= sqrt(rho)  # adjust for rho
            iso_x.append(p[0, 0] + np.pi)
            iso_y.append(p[1, 0])
        ax0.plot(iso_x, iso_y, lw=5, c=color, label="rho %.2f" % rho)
        ax1.plot(iso_x, iso_y, lw=5, c=color, label="rho %.2f" % rho)


    # contours are *point* based plots, so convert our bound into point
    # centers
    # cf = ax1.contourf(Qx[:-1, :-1] + Qdx/2.,
    #                   Qy[:-1, :-1] + Qdy/2., Qz2,
    #                   cmap=cmap)

    cf = ax1.imshow(Qz2, interpolation='bilinear', cmap=cmap,
                    origin='lower', extent=[x_range[0], x_range[1], y_range[0], y_range[1]],
                    vmax=0.0, vmin=Qz2.min(), aspect='auto')
    fig.colorbar(cf, ax=ax1)
    ax1.set_title(r"$\dot{V}(x,v)$")
    ax0.axis([x_range[0], x_range[1], y_range[0], y_range[1]])
    ax1.axis([x_range[0], x_range[1], y_range[0], y_range[1]])

    # adjust spacing between subplots so `ax1` title and `ax0` tick labels
    # don't overlap
    fig.tight_layout()

    plt.show()


class Dpoly(dict):

    def __add__(self, other):
        if isinstance(other, Dpoly):
            return Dpoly(poly_add(self, other))
        elif isinstance(other, float):
            ret = dict(self)
            tup = ret.keys()[0]
            tup_prime = tuple([0] * len(tup))
            if not ret.has_key(tup_prime):
                ret[tup_prime] = 0
            ret[tup_prime] += other
            return Dpoly(ret)
        raise NotImplementedError()

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, Dpoly):
            return Dpoly(poly_mult(self, other))
        elif isinstance(other, float):
            ret = dict(self)
            for key in self.keys():
                ret[key] = self[key] * other
            return Dpoly(ret)
        raise NotImplementedError()

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        ret = dict(self)
        for key in self.keys():
            ret[key] = -self[key]
        return Dpoly(ret)

    def __repr__(self):
        return string_from_dict(self)

    def __call__(self, *args):
        s=0.0
        for k in self.keys():
            s+=self[k] * reduce(
                lambda x,y: x*y, 
                [pow(args[i],j) for i, j in enumerate(list(k))]
                )
        return s


def test_poly_eval():
    dVp = Dpoly({(2,0):1.0, (1,1):3.0, (0,2):17.0})
    print "dVp", dVp
    print "dVp(0.,0.)", dVp(0.,0.)
    assert abs(dVp(0.,0.)-0.0)<1e-7
    print "dVp(1.,0.)", dVp(1.,0.)
    assert abs(dVp(1.,0.)-1.0)<1e-7
    print "dVp(2.,0.)", dVp(1.,0.)
    assert abs(dVp(2.,0.)-4.0)<1e-7
    print "dVp(0.,1.)", dVp(0.,1.)
    assert abs(dVp(0.,1.)-17.0)<1e-7
    print "dVp(1.,1.)", dVp(1.,1.)
    assert abs(dVp(1.,1.)-21.0)<1e-7







def play_with_poly_algebra(rho=0.5):
    print "rho is", rho
    V = Dpoly({(2,): 1.0})
    Vd = Dpoly({(2,): -1.0, (4,): 1.0})
    h = Dpoly({(2,): 1.0})
    print Vd
    print -Vd + rho
    print rho
    print rho + Vd
    print rho + V
    print rho - V
    print h * V
    print V * h
    print h * (rho + V)
    print h * (rho - V)
    print h * (V + rho)
    print h * (V - rho)
    print -Vd + 1.0
    print "test", -Vd
    print "test", -Vd + h * (V - rho)
    print "test", V
    print "test", h * V
    print "test", h * (V - rho)


def Htest(assert_feasible, rho=0.5):
    print "The following example tests the inclusion of an h polynomial."
    print "let Vd = -x^2 + x^4"
    print "let V = x^2"
    print "clearly Vd is negative until abs(x)>1"
    print "we want -Vd > 0"
    V = Dpoly({(2,): 1.0})
    Vd = Dpoly({(2,): -1.0, (4,): 1.0})
    h = Dpoly({(2,): 1.01})
    print "V =", V
    print "Vd =", Vd
    print "an early example estimate of h is", h
    print "rho is", rho
    print "poly we care about, h=1.01:", -Vd + h * (V - rho)

    h = Dpoly({(2, ): cvx.Variable(name="h_x2")})
    print h

    v_monomials = [(1, ), (2, )]
    h_monomials = [(1, )]

    Vmat = cvx.Semidef(len(v_monomials), name="Vmat")
    Hmat = cvx.Semidef(len(h_monomials), name="Hmat")

    obj = cvx.Minimize(cvx.norm(Vmat) + cvx.norm(Hmat))

    constraints = []

    poly_to_optimize = -Vd + h * (V - rho)

    print poly_to_optimize

    constraints.extend(setupSOSconstraint(
        Vmat, poly_to_optimize, v_monomials))
    constraints.extend(setupSOSconstraint(
        Hmat, h, h_monomials))
    prob = cvx.Problem(obj, constraints)
    prob.solve(solver=cvx.MOSEK, verbose=True)

    print prob.status
    assert (prob.status == assert_feasible)


def SOStest(V, Vd, assert_feasible, order=6, rho=0.5):
    V = Dpoly(V)
    Vd = Dpoly(Vd)
    print "V =", V
    print "Vd =", Vd
    print "rho is", rho
    h = Dpoly({
        (2, 0): cvx.Variable(name="h20"), 
        (1, 1): cvx.Variable(name="h11"), 
        (0, 2): cvx.Variable(name="h02"), 
        (4, 0): cvx.Variable(name="h40"), 
        (3, 1): cvx.Variable(name="h31")
        })
    print "poly we care about, h=1.01:", -Vd + h * (V - rho)

    if order==6:
        v_monomials = [(1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (3, 0)]
        h_monomials = [(1, 0), (0, 1), (2, 0), (2, 1)]
    if order==4:
        v_monomials = [(1, 0), (2, 0), (0, 1), (1, 1)]
        h_monomials = [(1, 0), (0, 1)]
    Vmat=cvx.Semidef(len(v_monomials), name = "Vmat")
    Hmat=cvx.Semidef(len(h_monomials), name = "Hmat")

    obj=cvx.Minimize(cvx.norm(Vmat) + cvx.norm(Hmat))

    constraints=[]

    poly_to_optimize=-Vd + h * (V - rho)

    print poly_to_optimize

    constraints.extend(setupSOSconstraint(
        Vmat, poly_to_optimize, v_monomials))
    constraints.extend(setupSOSconstraint(
        Hmat, h, h_monomials))
    prob=cvx.Problem(obj, constraints)
    try:
        prob.solve(solver = cvx.MOSEK, verbose = True)
        status=prob.status
    except cvx.error.SolverError:
        status="error"

    print status
    if status=="optimal":
        print Vmat.value
        print Hmat.value
    assert (status == assert_feasible)


def SOS4test(V, Vd, assert_feasible, rho=0.5):
    V = Dpoly(V)
    Vd = Dpoly(Vd)
    print "V =", V
    print "Vd =", Vd
    print "rho is", rho
    h = Dpoly({
        (2, 0): cvx.Variable(name="h20"), 
        (1, 1): cvx.Variable(name="h11"), 
        (0, 2): cvx.Variable(name="h02")
        })
    print "poly we care about:", -Vd + h * (V - rho)


    v_monomials = [(1, 0), (0, 1), (2, 0), (1, 1), (0,2)]
    h_monomials = [(1, 0), (0, 1)]
    Vmat=cvx.Semidef(len(v_monomials), name = "Vmat")
    Hmat=cvx.Semidef(len(h_monomials), name = "Hmat")

    obj=cvx.Minimize(cvx.norm(Vmat) + cvx.norm(Hmat))

    constraints=[]

    poly_to_optimize=-Vd + h * (V - rho)

    print poly_to_optimize

    constraints.extend(setupSOSconstraint(
        Vmat, poly_to_optimize, v_monomials))
    constraints.extend(setupSOSconstraint(
        Hmat, h, h_monomials))
    prob=cvx.Problem(obj, constraints)
    try:
        prob.solve(solver = cvx.MOSEK, verbose = True)
        status=prob.status
    except cvx.error.SolverError:
        status="error"

    print status
    if status=="optimal":
        print Vmat.value
        print np.linalg.eigh(Vmat.value)
        print Hmat.value
        print np.linalg.eigh(Hmat.value)

    assert (status == assert_feasible)

def SOS4check(V, Vd, rho=0.5):
    V = Dpoly(V)
    Vd = Dpoly(Vd)
    h = Dpoly({
        (2, 0): cvx.Variable(name="h20"), 
        (1, 1): cvx.Variable(name="h11"), 
        (0, 2): cvx.Variable(name="h02")
        })
    v_monomials = [(1, 0), (0, 1), (2, 0), (1, 1), (0,2)]
    h_monomials = [(1, 0), (0, 1)]
    Vmat=cvx.Semidef(len(v_monomials), name = "Vmat")
    Hmat=cvx.Semidef(len(h_monomials), name = "Hmat")

    obj=cvx.Minimize(cvx.norm(Vmat) + cvx.norm(Hmat))

    constraints=[]

    poly_to_optimize=-Vd + h * (V - rho)


    constraints.extend(setupSOSconstraint(
        Vmat, poly_to_optimize, v_monomials))
    constraints.extend(setupSOSconstraint(
        Hmat, h, h_monomials))
    prob=cvx.Problem(obj, constraints)
    try:
        prob.solve(solver = cvx.MOSEK, verbose = False)
        status=prob.status
    except cvx.error.SolverError:
        status="error"

    print status
    if status=="optimal":
        return 1.0
    if status=="optimal_inaccurate":
        return 0.75
    if status=="infeasible":
        return 0.0
    if status=="error":
        return 0.0
    assert (False)


def test_order_six():
    plotrho([100.,120.], color = 'r' if prob.status != "optimal" else 'g', order=6)
    Vpoly=nosympy(listtodict(sparse_taylor(Vprime, dX, dX0, N=2, tol=1e-9)))
    Vdotpoly=nosympy(listtodict(sparse_taylor(dV, dX, dX0, N=6, tol=1e-9)))
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 0.125)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 0.25)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 0.5)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 1.0)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 2.5)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 5.0)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 7.5)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 10.0)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 15.0)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 30.0)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 60.0)
    SOStest(Vpoly, Vdotpoly, "optimal", order=6,rho = 100.0)
    SOStest(Vpoly, Vdotpoly, "optimal_inaccurate", order=6,rho = 120.0)
    SOStest(Vpoly, Vdotpoly, "optimal_inaccurate", order=6,rho = 140.0)
    SOStest(Vpoly, Vdotpoly, "optimal_inaccurate", order=6,rho = 160.0)
    SOStest(Vpoly, Vdotpoly, "optimal_inaccurate", order=6,rho = 180.0)
    SOStest(Vpoly, Vdotpoly, "optimal_inaccurate", order=6,rho = 200.0)
    SOStest(Vpoly, Vdotpoly, "optimal_inaccurate", order=6,rho = 220.0)
    SOStest(Vpoly, Vdotpoly, "error", order=6,rho = 235.0)
    SOStest(Vpoly, Vdotpoly, "error", order=6,rho = 260.0)
    SOStest(Vpoly, Vdotpoly, "error", order=6,rho = 300.0)
    SOStest(Vpoly, Vdotpoly, "infeasible", order=6,rho = 340.0)
    SOStest(Vpoly, Vdotpoly, "infeasible", order=6,rho = 360.0)
    SOStest(Vpoly, Vdotpoly, "infeasible", order=6,rho = 380.0)
    SOStest(Vpoly, Vdotpoly, "infeasible", order=6,rho = 400.0)
    SOStest(Vpoly, Vdotpoly, "infeasible", order=6,rho = 450.0)
    SOStest(Vpoly, Vdotpoly, "infeasible", order=6,rho = 490.0)
    SOStest(Vpoly, Vdotpoly, "infeasible", order=6,rho = 1000.0)

NsectN=15
def Nsect(rlow, rhigh, Vpoly, Vdotpoly):
    rhos = np.linspace(rlow,rhigh,NsectN)
    results = [1.0]+[SOS4check(Vpoly, Vdotpoly, rho=r) for r in rhos[1:-1]]+[0.0]
    index_of_first_unsure_result=NsectN
    for i in range(NsectN):
        if results[i]!=1.0:
            index_of_first_unsure_result=i-1
            break
    index_of_last_unsure_result=0.0
    for i in reversed(range(NsectN)):
        if results[i]!=0.0:
            index_of_last_unsure_result=i+1
            break
    print rhos
    print results
    print "index_of_last_unsure_result", index_of_last_unsure_result
    print "index_of_first_unsure_result", index_of_first_unsure_result
    assert index_of_last_unsure_result>index_of_first_unsure_result

    new_rlow = rhos[index_of_first_unsure_result]
    new_rhigh = rhos[index_of_last_unsure_result]
    if new_rlow==rlow and new_rhigh==rhigh:
        return rlow, rhigh
    print "next Nsect(%.2f, %.2f)"%(new_rlow,new_rhigh)
    return Nsect(new_rlow,new_rhigh, Vpoly, Vdotpoly)




def test_order_four():
    Vpoly=nosympy(listtodict(sparse_taylor(Vprime, dX, dX0, N=2, tol=1e-9)))
    Vdotpoly=nosympy(listtodict(sparse_taylor(dV, dX, dX0, N=4, tol=1e-9)))
    SOS4test(Vpoly, Vdotpoly, "optimal", rho = 10.0)
    SOS4test(Vpoly, Vdotpoly, "optimal", rho = 15.0)
    SOS4test(Vpoly, Vdotpoly, "optimal", rho = 20.0)
    SOS4test(Vpoly, Vdotpoly, "optimal_inaccurate", rho = 25.0)
    SOS4test(Vpoly, Vdotpoly, "infeasible", rho = 30.0)
    SOS4test(Vpoly, Vdotpoly, "infeasible", rho = 35.0)
    SOS4test(Vpoly, Vdotpoly, "infeasible", rho = 40.0)

    SOS4test(Vpoly, Vdotpoly, "infeasible", rho = 100.0)
    # SOS4test(Vpoly, Vdotpoly, "optimal", rho = 43.0)
    SOS4test(Vpoly, Vdotpoly, "infeasible", rho = 300.0)


    plotrho([15.,30.], color = 'g', order=4)











if __name__ == '__main__':

    PrintUTHeader()
    test_poly_eval()

    # some notes about using sympy:
    # sp.N is a function for numerically approximating results as floats
    # sp.lambda converts an expression into a function, Xdot[1][0] is a
    # function of x

    print "here we practice using sympy's taylor, diff, and my new part_taylor function"
    print "taylor expansion of f[1][0] w.r.t. x:", [sp.N(res) for res in taylor(sp.Lambda(x, f[1][0]), pi, 6)]
    print "clearly the taylor expansion function has serious problems with the included v term."
    # now looking into direct differentiation

    # Generate linear approximation
    X0=np.array([[pi], [0.0]])  # top of arc
    u0=np.array([[0.0]])  # no current required at this point

    print "using diff to find the fifth term independantly:", f[1][0].diff(x, 5).subs(x, pi) / factorial(5)
    print "using diff to find the sixth term independantly:", f[1][0].diff(x, 6).subs(x, pi) / factorial(6)

    A, B=lin(f, X, X0)

    print
    print "Linearized system matrix A:\n", A
    print "Linearized input matrix B:\n", B
    print

    print
    print "We also require that the nominal input and state are actually an equilibrium"
    print "0 = f(x0)+g(x0)u0"
    print(f + g * u)
    print " = "
    res=np.array([
        [sp.N((f + g * u)[0, 0].subs(x, X0[0, 0]).subs(v, X0[1, 0]).subs(u, u0))],
        [sp.N((f + g * u)[1, 0].subs(x, X0[0, 0]).subs(v, X0[1, 0]).subs(u, u0))]])
    if sum(abs(res)) < 1e-7:
        print np.array([[0], [0]])
        print "close enough, anyway"
    else:
        print res
    print

    print
    print "Now we define the cost matrices for LQR"
    print "A:\n", A
    print "B:\n", B
    print "Q:\n", Q
    print "R:\n", R
    print
    # use slycot, wrapped inside python controls, to compute the algebraic
    # riccati equation's solution
    print "\nDo LQR. the control module is old, and generates warnings:"

    S, L, K=ctrl.care(A, B, Q, R)
    print "check: A'S+SA-SBRinvB'S+Q=", A.T.dot(S) + S.dot(A) - S.dot(B).dot(np.linalg.inv(R)).dot(B.T).dot(S) + Q
    # exit()

    print
    print "Symmetric, positive definite, quadratic-form-lypunov-function matrix S:\n", S
    print "Eigenvalues of the new closed loop", L
    print "Full state feedback matrix K", K
    print "Eigenvalues of S", np.linalg.eigh(S)[0]
    print "Eigenvalues of Z", np.linalg.eigh(S.dot(B).dot(np.linalg.inv(R)).dot(B.T).dot(S) + Q)[0]
    print

    print
    print "now we need to setup the feasibility problem"
    rho=80.0
    print "assume rho is", rho
    dx, dv, du=sp.Symbol("dx"), sp.Symbol("dv"), sp.Symbol("du")
    dX=np.array([[dx], [dv]])
    dX0=np.array([[0], [0]])

    V=((X - X0).T.dot(S).dot((X - X0)))[0, 0]
    Vprime=V.subs(x, X0[0, 0] + dx).subs(v, X0[1, 0] + dv)
    print "Lyapunov function V(x,v) =", V
    print "linearizing about the operating point,"
    Vpoly=sparse_taylor(Vprime, dX, dX0, N = 9, tol = 1e-9)
    print "V ~=", string_from_taylor(sparse_taylor(Vprime, dX, dX0, N=9, tol=1e-9))
    print "raw poly repr is V=", Vpoly
    print "technically, this has just confirmed that we can represent lyapunov functions in polynomial format."
    print

    print
    print "Now we look at the derivative of the lyapunov function, as a taylor approximate"
    print "Derivative is", Vprime.diff(dx) * sp.Symbol("dxdt") + Vprime.diff(dv) * sp.Symbol("dvdt")
    print "where dxdt =", (f[0, 0] + g[0, 0] * u).subs(x, X0[0, 0] + dx).subs(v, X0[1, 0] + dv)
    print "  and dvdt =", (f[1, 0] + g[1, 0] * u).subs(x, X0[0, 0] + dx).subs(v, X0[1, 0] + dv)
    dV=V.diff(x) * (f[0, 0] + g[0, 0] * u) + \
        V.diff(v) * (f[1, 0] + g[1, 0] * u)
    print "thus dot{V} =", dV.subs(x, X0[0, 0] + dx).subs(v, X0[1, 0] + dv)
    dV=dV.subs(u, (u0 - (K.dot(dX)))
                 [0, 0]).subs(x, X0[0, 0] + dx).subs(v, X0[1, 0] + dv)
    print "  or dot{V} =", dV
    print "which we can approximate as a big polynomial in dx and dv"
    print "dot{V} ~=", string_from_taylor(sparse_taylor(dV, dX, dX0, N=4, tol=1e-9))
    print

    print
    print "our first exercize in sum of squares verification might be to show that the second order approximation"

    print "dov{V} ~=", string_from_taylor(sparse_taylor(dV, dX, dX0, N=2, tol=1e-9))
    print "is sum of squares. Actually, we already know it is since it satisifes the riccati equation."
    print "A'S + S A - SBR^{-1}B'S +Q = 0"
    print "starting with A'S + S A < 0"
    print "A <- A - BK"
    print "A'S - K'B'S + SA - SBK < 0"
    print "K = R^{-1}B'S"
    print "A'S + SA - 2 * SBR^{-1}B'S"
    print "substituting the riccati equation"
    print "0 <  SBR^{-1}B'S + Q "
    Z=S.dot(B).dot(np.linalg.inv(R)).dot(B.T).dot(S) + Q
    print "0 < ...\n", Z
    print "where this matrix, call it Z, is such that dot{V} = -[dx, dv]*Z*[dx; dv]"
    print "which we can verify using sympy.\ndot{V} =", -dX.T.dot(Z).dot(dX)[0, 0]
    print "dov{V} ~=", string_from_taylor(sparse_taylor(dV, dX, dX0, N=2, tol=1e-9), precision = 5)
    print

    print
    print "Now let us consder a polynomial I came up with, in one variable."
    mypoly={(2,): 1.0, (3,): 0.0, (4,): 2.000, (5,): 1.0, (6,): 1.0}
    print string_from_dict(mypoly)
    print "This can be expressed using the monomial vector"
    monomials=[(1,), (2,), (3,)]
    print[string_from_dict({k: 1.0}) for k in monomials]

    P=S

    S=cvx.Semidef(3)
    obj=cvx.Minimize(cvx.norm(S))
    dct={}
    for i in range(3):
        for j in range(3):
            monom=tuple([monomials[i][q] + monomials[j][q]
                           for q in range(1)])
            if not dct.has_key(monom):
                dct[monom]=[]
            dct[monom].append((i, j))
    constraints=[]
    print "here are the 'poly constraints', as I like to call them:"
    for monom in dct.keys():
        if mypoly.has_key(monom):
            coeff=mypoly[monom]
        else:
            coeff=0.0
        constraints.append(coeff == sum([S[i, j] for i, j in dct[monom]]))
        print "\t", coeff, "=", "+".join(["S[%d,%d]" % (i, j) for i, j in dct[monom]])
    prob=cvx.Problem(obj, constraints)
    prob.solve(solver = cvx.MOSEK)
    assert(prob.status == "optimal")
    # prob.solve(solver=cvx.CVXOPT)
    # prob.solve(solver=cvx.SCS)

    print S.value
    print np.linalg.eigh(S.value)
    print

    print
    print "now we scale up the difficulty: coefficients which are cvx variables."
    mypoly={(2,): cvx.Variable(), (3,): 0.0, (4,)
             : cvx.Variable(), (5,): 1.0, (6,): cvx.Variable()}
    print string_from_dict(mypoly)
    print "This can be expressed using the monomial vector"
    monomials=[(1,), (2,), (3,)]
    print[string_from_dict({k: 1.0}) for k in monomials]

    S=cvx.Semidef(3)
    obj=cvx.Minimize(cvx.norm(S))
    dct={}
    for i in range(3):
        for j in range(3):
            monom=tuple([monomials[i][q] + monomials[j][q]
                           for q in range(1)])
            if not dct.has_key(monom):
                dct[monom]=[]
            dct[monom].append((i, j))
    constraints=[]
    print "here are the 'poly constraints', as I like to call them:"
    for monom in dct.keys():
        if mypoly.has_key(monom):
            coeff=mypoly[monom]
        else:
            coeff=0.0
        constraints.append(coeff == sum([S[i, j] for i, j in dct[monom]]))
        print "\t", coeff, "=", "+".join(["S[%d,%d]" % (i, j) for i, j in dct[monom]])
    constraints.append(mypoly[(2,)] == 1.0)
    constraints.append(mypoly[(4,)] == 1.5)
    constraints.append(mypoly[(6,)] == 1.0)
    prob=cvx.Problem(obj, constraints)
    prob.solve(solver = cvx.MOSEK)
    assert(prob.status == "optimal")
    # prob.solve(solver=cvx.CVXOPT)
    # prob.solve(solver=cvx.SCS)

    print S.value
    print


    def Ptest(a, b, assert_feasible):
        print "The following example tests the inclusion of an h polynomial."
        print "let Vd = -x^2 + x^4"
        print "let V = x^2"
        print "clearly Vd is negative until abs(x)>1"
        print "we want -Vd > 0"
        print "but we specify to -Vd + h(x)(V-0.5) > 0, h(x)>0"
        print "if h(x) = x^2, then "
        print "-x^4+x^2+x^4-0.5x^2 > 0"
        print "clearly illustrates that Vd is neagative for V<0.5"


        Vdotpoly= {(2,): a, (4,): b}
        v_monomials=[(1, ), (2, )]

        Vmat=cvx.Semidef(len(v_monomials), name = "Vmat")
        obj=cvx.Minimize(cvx.norm(Vmat))
        constraints=[]
        constraints.extend(setupSOSconstraint(
            Vmat, Vdotpoly, v_monomials))

        prob=cvx.Problem(obj, constraints)
        prob.solve(solver = cvx.MOSEK, verbose = True)
        print prob.status
        assert (prob.status == assert_feasible)




    Ptest(1.0, 1.0, 'optimal')
    Ptest(1.0, 100.0, 'optimal')
    Ptest(100.0, 1.0, 'optimal')
    Ptest(1000.0, 1000.0, 'optimal')
    Ptest(-1000.0, 1000.0, 'infeasible')
    Ptest(-1.0, 1.0, 'infeasible')
    Ptest(-1.0e-2, 1.0, 'infeasible')
    Ptest(-1.0e-4, 1.0, 'infeasible')
    Ptest(-1.0e-4, 1.0e4, 'infeasible')  # only works with MOSEK
    Htest('infeasible', rho = 2.0)
    bad_values_violate_assertion=False
    try:
        Htest('infeasible', rho = 0.5)
    except AssertionError:
        bad_values_violate_assertion=True
    assert bad_values_violate_assertion
    Htest('optimal', rho = 0.5)
    Htest('optimal', rho = 0.1)
    Htest('optimal', rho = 0.9)
    Htest('infeasible', rho = 1.1)


    print

    print


    Vpoly=nosympy(listtodict(sparse_taylor(Vprime, dX, dX0, N=2, tol=1e-9)))
    Vdotpoly=nosympy(listtodict(sparse_taylor(dV, dX, dX0, N=4, tol=1e-9)))
    rholow, rhohigh = Nsect(1.0, 100., Vpoly, Vdotpoly)
    plotrho([rholow,rhohigh], color = 'g', order=4)
    # test_order_four()
    # test_order_six()

    ##########################################################################

    exit()

    ##########################################################################

    # exit()
    print

    print "now we go to checking a particular rho,", rho
    print "we check that -dot{V}+h(x,v)*(V(x,v)-rho) is SOS"
    print "and that h(x,v) is SOS"
    Vpoly=listtodict(sparse_taylor(Vprime, dX, dX0, N=2, tol=1e-9))
    print "V ~=", string_from_dict(Vpoly)
    Vdotpoly=listtodict(sparse_taylor(dV, dX, dX0, N=4, tol=1e-9))
    print "dot{V} ~=", string_from_dict(Vdotpoly)


    h={(2, 0): cvx.Variable(name = "h1"), (1, 1): cvx.Variable(
        name = "h2"), (0, 2): cvx.Variable(name = "h3")}
    print "h =", string_from_dict(h)
    print Vpoly
    print Vdotpoly
    print h

    print "Total expression: -Vdot + h*(V-rho) is epsilon-SOS and h is epsilon-SOS"
    VasPoly=nosympy(dict(Vpoly))
    print "V as a polynomial"

    VlessRho=dict(VasPoly)
    VlessRho[0, 0]=-rho
    print "V-rho as a polynomial"

    print VlessRho
    htimesVlessRho=poly_mult(VlessRho, h)
    expr1=poly_add(poly_add(poly_mult(nosympy(Vdotpoly), {
                     (0, 0): -1.0}), htimesVlessRho), {(2, 0): -1e-3, (0, 2): -1e-3})
    print "h(V-rho) has coefficients for the following monomial orders:", htimesVlessRho.keys()
    print "-dot{V}+h(V-rho) has the following monomial orders:", expr1.keys()
    print "for this we use the monomial vector:"
    v_monomials=[(1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
    h_monomials=[(1, 0), (0, 1)]
    print[string_from_dict({k: 1.0}, names = ['dx', 'dv']) for k in v_monomials]
    nummon=len(v_monomials)

    Vmat=cvx.Semidef(len(v_monomials), name = "Vmat")
    Hmat=cvx.Semidef(len(h_monomials), name = "Vmat")
    obj=cvx.Minimize(cvx.norm(Vmat))
    constraints=[]
    constraints.extend(setupSOSconstraint(
        Vmat + 0.0 * np.eye(len(v_monomials)), htimesVlessRho, v_monomials))
    constraints.extend(setupSOSconstraint(
        Hmat + 0.0 * np.eye(len(h_monomials)), h, h_monomials))
    prob=cvx.Problem(obj, constraints)
    prob.solve(solver = cvx.SCS, verbose = True)
    print Vmat.value
    print Hmat.value
    if Hmat.value != None:
        print np.linalg.eigh(Hmat.value + 1e-4 * np.eye(len(h_monomials)))
    print prob.status
    plotrho(rho, color = 'r' if prob.status != "optimal" else 'g')
    # assert(prob.status=="optimal")
    exit()
