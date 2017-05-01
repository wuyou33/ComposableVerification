% Define variables


P1=[ 100,0;0,1];
P2= [ 6.0000   -0.5000
   -0.5000    5.5000];
U1M=[0;1];
U2M=[.5257;.8507];

U1m=[1/100;0];
U2m=[ -0.8507  ;  0.5257];
  
    
% % % looking for the touching hyperplane
x1=sdpvar(2,1);
x2 = sdpvar(2,1);
% variable x2(2);



% Define constraints and objective

Constraints = [[1,x1';x1,inv(P1)]>=1e-7*eye(3),
[1,x2';x2,inv(P2)]>=1e-7*eye(3),P1*x1==P2*x2];
Objective = -(U1M'*x1+U2M'*x2);
% Set some options for YALMIP and solver
% options = sdpsettings('verbose',1,'solver','quadprog','quadprog.maxiter',100);
% Solve the problem
sol = optimize(Constraints,Objective);
% Analyze error flags
if sol.problem == 0
 % Extract and display value
 x1 = value(x1)
 x2 = value(x2)
else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem)
end