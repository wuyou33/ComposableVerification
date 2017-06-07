degV=4;

x=msspoly('x',16);
y=msspoly('y',3);


x1=x(1);
x2=x(2);
x3=x(3);
x4=x(4);
x5=x(5);
x6=x(6);
x7=x(7);
x8=x(8);
x9=x(9);
x10=x(10);
x11=x(11);
x12=x(12);
x13=x(13);
x14=x(14);
x15=x(15);
x16=x(16);

X1 = [x1,x3,x7,x8,x9]';
X2 = [x2,x6,x10,x15,x16]';
X3 = [x4,x5,x11,x12,x13,x14]';


dx1dt = -1.0*x1^2-0.167*x1*x3-1.16*x1-0.193*x3;
dx3dt = -1.0*x3^2-0.397e-1*x3;
dx7dt = -0.143*x1*x7-1.0*x7^2-0.388*x8*x7-0.163*x1-1.14*x7-0.443*x8;
dx8dt = -1.0*x8^2-0.198*x8;
dx9dt = -0.146*x9*x3-1.0*x9^2-0.107*x3-0.730*x9; 
dX1dt=[dx1dt;dx3dt;dx7dt;dx8dt;dx9dt];
g9_12=-0.994e-2*x12-0.136e-1*x9*x12;


dx2dt = -1.0*x2^2-0.113*x10*x2-1.26*x2-0.142*x10;
dx6dt = -1.0*x6^2-0.313*x6; 
dx10dt = -0.419*x6*x10-1.0*x10^2-0.822*x6-1.96*x10;
dx15dt = -0.796*x15*x6-0.230*x15*x10-1.0*x15^2-1.36*x6-0.394*x10-1.71*x15;
dx16dt = -0.527*x15*x16-1.0*x16^2-0.356*x15-0.676*x16;
dX2dt = [dx2dt,dx6dt,dx10dt,dx15dt,dx16dt]';
g2_4_5=-1.22*x4-0.333*x5-0.969*x4*x2-0.264*x5*x2;
g46=-0.866*x6*x4-0.271*x4;
g15_3=-0.548*x15*x3-0.936*x3;
g16_5=-0.407*x16*x5-0.275*x5;


dx4dt = -1.0*x4^2-0.926e-1*x4;
dx5dt = -0.847*x4*x5-1.0*x5^2-1.77*x4-2.09*x5;
dx11dt = -1.0*x11^2-0.391e-1*x11;
dx12dt = -1.0*x12^2-0.244*x12;
dx13dt = -0.211*x5*x13-1.0*x13^2-0.193*x5-0.915*x13;
dx14dt = -0.709*x14*x12-1.0*x14^2-1.09*x12-1.54*x14;
dX3dt = [dx4dt,dx5dt,dx11dt,dx12dt,dx13dt,dx14dt]';
g11_10=-0.947e-1*x11*x10-0.370e-2*x10;
g12_1=-0.801*x1*x12-0.195*x1;
g14_9=-0.571*x14*x9-0.881*x9;


% construct V
V1Xmonom = monomials(X1,0:degV);
V2Xmonom = monomials(X2,0:degV);
V3Xmonom = monomials(X3,0:degV);


% partition 1

prog = spotsosprog;
prog = prog.withIndeterminate([X1;y(1)]);

%V condition
[prog,V1] = prog.newFreePoly(V1Xmonom);
prog = prog.withSOS(V1);
PV1PX1 = diff(V1,X1);
V1dot = diff(V1,X1)*dX1dt;
V1dotcomp=-V1dot+2*X1'*PV1PX1'*(0.136e-1*x9)*y(1)+y(1)^2-(-0.936*x3)^2-(-0.195*x1)^2-(x1)^2-(-0.881*x9)^2-(x9)^2;
prog = prog.withSOS(-V1dot);
prog = prog.withSOS(-V1dotcomp);

options = spot_sdp_default_options();
sol = prog.minimize(0,@spot_mosek,options);

if ~sol.isPrimalFeasible
  error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
  error('Problem looks dual infeasible. It is probably unbounded. ');
end

V1 = sol.eval(V1);
partition2=0;
partition3=0;
% partition 2
if (partition2==1)
	prog = spotsosprog;
	prog = prog.withIndeterminate([X2;y(2)]);
	[prog,V2] = prog.newFreePoly(V2Xmonom);
	prog = prog.withSOS(V2);
	PV2PX2 = diff(V2,X2);
	V2dot = diff(V2,X2)*dX2dt;
	prog = prog.withSOS(-(V2dot+2*X2'*PV2PX2'*y(2)-y(2)'*y(2))); 
	sol = prog.minimize(0,@spot_mosek,options);

	if ~sol.isPrimalFeasible
	  error('Problem looks primal infeasible');
	end

	if ~sol.isDualFeasible
	  error('Problem looks dual infeasible. It is probably unbounded. ');
	end
	V2= sol.eval(V2);
end

if (partition3==1)
% partition 3
prog = spotsosprog;
prog = prog.withIndeterminate([X3;y(3)]);
[prog,V3] = prog.newFreePoly(V3Xmonom);
prog = prog.withSOS(V3);
PV3PX3 = diff(V3,X3);
V3dot = diff(V3,X3)*dX3dt;
prog = prog.withSOS(-(V3dot+2*X3'*PV3PX3'*y(3)-y(3)'*y(3))); 
sol = prog.minimize(0,@spot_mosek,options);

if ~sol.isPrimalFeasible
  error('Problem looks primal infeasible');
end

if ~sol.isDualFeasible
  error('Problem looks dual infeasible. It is probably unbounded. ');
end
V3= sol.eval(V3);
end

disp(V1);
disp(V2);
disp(V3)

% rho = double(sol.eval(rho));