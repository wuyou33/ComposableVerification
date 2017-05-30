degV1=2;
degV2=2;

x12=msspoly('x12',2);
x3=msspoly('x3',1);
x1=x12[1];
x2=x12[2];

V1monom = monomials(x12,0:degV1);
V2monom = monomials(x3,0:degV2);

prog.withSOS()
prog = spotsosprog;
prog = prog.withIndeterminate(x12);
dx1dt = x1*(-4*x1 + 5*x2 +1);
dx2dt = x2*(-8*x1 -12*x2 - 8*x3 + 6);
dx3dt = x3*(5*x2 - 4*x3 + 1);

f12=[dx1dt;dx2dt];
f3=dx3dt;

V1=prog.withSOS()

  % construct Vdot
Vdot = clean(diff(V,x)*f);
  
  % construct slack var
  [prog,sigma1] = prog.newPos(1);

  % setup SOS constraints
  prog = prog.withSOS(-Vdot + L1*(V - 1) - sigma1*V);
  prog = prog.withSOS(L1);

  % run SeDuMi/MOSEK and check output
  solver = options.solver;
  options = spot_sdp_default_options();
  sol = prog.minimize(-sigma1,solver,options);
  
  if sol.status == spotsolstatus.STATUS_SOLVER_ERROR
    error('The solver threw an internal error.');
  end
   if ~sol.isPrimalFeasible
      error('Problem looks primal infeasible');
  end
  
  if ~sol.isDualFeasible
      error('Problem looks dual infeasible. It is probably unbounded. ');
  end

  L1 = sol.eval(L1);
  sigma1 = sol.eval(sigma1);
end