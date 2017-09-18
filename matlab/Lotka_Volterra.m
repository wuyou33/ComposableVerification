% the dynamics
% dx1/dt = -1.0*x1^2-0.167*x1*x3-1.16*x1-0.193*x3;
% dx2/dt = -1.0*x2^2-0.969*x4*x2-0.264*x5*x2-0.113*x10*x2-1.26*x2-1.22*x4-0.333*x5-0.142*x10;
% dx3/dt = -1.0*x3^2-0.397e-1*x3;
% dx4/dt = -1.0*x4^2-0.926e-1*x4;
% dx5/dt = -0.847*x4*x5-1.0*x5^2-1.77*x4-2.09*x5;
% dx6/dt = -0.866*x6*x4-1.0*x6^2-0.271*x4-0.313*x6;
% dx7/dt = -0.143*x1*x7-1.0*x7^2-0.388*x8*x7-0.163*x1-1.14*x7-0.443*x8;
% dx8/dt = -1.0*x8^2-0.198*x8;
% dx9/dt = -0.146*x9*x3-1.0*x9^2-0.136e-1*x9*x12-0.107*x3-0.730*x9-0.994e-2*x12;
% dx10/dt = -0.419*x6*x10-1.0*x10^2-0.822*x6-1.96*x10;
% dx11/dt = -0.947e-1*x11*x10-1.0*x11^2-0.370e-2*x10-0.391e-1*x11;
% dx12/dt = -0.801*x1*x12-1.0*x12^2-0.195*x1-0.244*x12;
% dx13/dt = -0.211*x5*x13-1.0*x13^2-0.193*x5-0.915*x13;
% dx14/dt = -0.571*x14*x9-0.709*x14*x12-1.0*x14^2-0.881*x9-1.09*x12-1.54*x14;
% dx15/dt = -0.548*x15*x3-0.796*x15*x6-0.230*x15*x10-1.0*x15^2-0.936*x3-1.36*x6-0.394*x10-1.71*x15;
% dx16/dt = -0.407*x16*x5-0.527*x15*x16-1.0*x16^2-0.275*x5-0.356*x15-0.676*x16


function [time1,time2,time3,time] = Lotka_Volterra()
checkDependency('spotless');
checkDependency('mosek');

x=msspoly('x',16);
y=msspoly('y',6);
z=msspoly('z',5);

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

V1=0.559*x8^2+0.343*x8*x7+0.711*x7^2+0.435e-1*x8*x1+0.619e-2*x1*x7+0.725*x1^2+.222e-1*x8*x3-0.239e-1*x7
         0.238*x1*x3+.153*x3^2+0.741e-1*x9*x8+0.331e-1*x9*x7+0.774e-1*x9*x1+0.214*x9*x3+0.755*x9^2;

V2=0.546*x15^2+0.114e-1*x15*x16+0.789*x16^2+0.388*x15*x6-0.151e-1*x16*x6+1.00*x6^2-0.101*x15*x10+0.372e
         0.228*x6*x10+0.611*x10^2+0.175e-1*x15*x2+0.215e-1*x16*x2+0.138e-1*x6*x2-0.344e-2*x10*x2+0.795*x2^2;

V3=0.719*x4^2+0.500*x4*x5+0.338*x5^2+0.238e-2*x14*x4-0.443e-2*x14*x5+0.533*x14^2+0.653e-1*x12*x4+
                  0.126e-1*x12*x5+0.438*x14*x12+0.755*x12^2-0.139*x13*x4+0.163e-1*x5*x13+0.255e-1*x14*x13+0.305e-1*x13
                  0.784*x13^2+0.220e-1*x11*x4+0.297e-3*x11*x5-0.639e-3*x14*x11+0.761e-2*x12*x11-0.114e-3*x13*x11+0.806e-1*x11^2
% the partition
X1 = [x1,x3,x7,x8,x9]';
X2 = [x2,x6,x10,x15,x16]';
X3 = [x4,x5,x11,x12,x13,x14]';
epsi=1e-7;

f1=[-1.0*x1^2-0.167*x1*x3-1.16*x1-0.193*x3;
-1.0*x3^2-0.397e-1*x3;
-0.143*x1*x7-1.0*x7^2-0.388*x8*x7-0.163*x1-1.14*x7-0.443*x8;
-1.0*x8^2-0.198*x8;
-0.146*x9*x3-0.107*x3-0.730*x9-1.0*x9^2;];
g13=[zeros(4,1);-0.994e-2-0.136e-1*x9]; %(5,1)
h13=[x12]; %(1,1)
g12=0; %1,1
h12=0; %1,1

f2=[-1.0*x2^2-0.113*x10*x2-1.26*x2-0.142*x10;
-1.0*x6^2-0.313*x6;
 -0.419*x6*x10-1.0*x10^2-0.822*x6-1.96*x10;
-0.796*x15*x6-0.230*x15*x10-1.0*x15^2-1.36*x6-0.394*x10-1.71*x15;
-0.527*x15*x16-1.0*x16^2-0.356*x15-0.676*x16];
g21=[0;0;0;-0.548*x15-0.936;0];%5,1
h21=[x3];%1,1
g23=[x2,1,0,0,0,0;
	0,0,-.271,0,0,-0.866*x6;
	0,0,0,0,0,0
	0,0,0,0,0,0;
	0,0,0,0,-0.407*x16-0.275,0];% (5,6)
h23=[-0.969*x4-0.264*x5;-1.22*x4-0.333*x5;x4;0;x5;x4];% (6,1)


f3=[ -1.0*x4^2-0.926e-1*x4;
-0.847*x4*x5-1.0*x5^2-1.77*x4-2.09*x5;
-1.0*x11^2-0.391e-1*x11;
-1.0*x12^2-0.244*x12;
-0.211*x5*x13-1.0*x13^2-0.193*x5-0.915*x13;
-0.709*x14*x12-1.0*x14^2-1.09*x12-1.54*x14;];
g31=[0,0;0,0;0,0;-0.801*x12-0.195,0;0,0;0,-0.571*x14-0.881];  %(6,2)
h31=[x1;x9]; %(2,1)
g32=[0;0;-0.947e-1*x11-0.370e-2;0;0;0]; %(6,1)
h32=[x10]; %(1,1)

time=0;

% finding V1
tic
prog = spotsosprog;
prog = prog.withIndeterminate([X1;z]);
Vmonom = monomials(X1,0:2);
[prog,V] = prog.newSOSPoly(Vmonom);
[prog,slacks] = prog.newPos(2);
q1=(diff(V,X1)*f1)-.5*(h21'*h21+h31'*h31);
prog=prog.withSOS(-q1-slacks(1));
s1=diff(V,X1)*f1*g13;
prog=prog.withSOS(-q1+2*z'*s1+z'*z-slacks(2));
options = spot_sdp_default_options();
options.verbose=0;
sol=prog.minimize(-sum(slacks),@spot_mosek,options);
V1=sol.eval(V)
slacks=sol.eval(slacks)
time1=toc;

% finding V2
tic
prog = spotsosprog;
prog = prog.withIndeterminate([X2;y]);
Vmonom = monomials(X2,0:2);
[prog,V] = prog.newFreePoly(Vmonom);
prog=prog.withSOS(V);

q2=clean(diff(V,X2)*f2)-(h12'*h12+h32'*h32)/2;
prog=prog.withSOS(-q2);

s2=.5*(diff(V,X2)*f2*g21);
prog=prog.withSOS(-q2+2*y*s2+y^2);
options = spot_sdp_default_options();
options.verbose=0;
sol=prog.minimize(-0,@spot_mosek,options);
V2=sol.eval(V);
time2=toc;

% % % fidning V3
% tic
% prog = spotsosprog;
% prog = prog.withIndeterminate([X3;y]);
% Vmonom = monomials(X3,0:2);
% [prog,V] = prog.newFreePoly(Vmonom);
% prog=prog.withSOS(V-epsi);

% q3=clean(diff(V,X3)*f3)-(h13'*h13+h23'*h23)/2;
% prog=prog.withSOS(-q3+epsi);

% s3=.5*(diff(V,X3)*f3*g32);
% prog=prog.withSOS(-q3+2*y*s3+y^2);

% options = spot_sdp_default_options();
% options.verbose=0;
% sol=prog.minimize(-0,@spot_mosek,options);
% V3=sol.eval(V);
% time3=toc;

% time=time1+time2+time3;
end


% the dynamics grouped by partition
% [-1.0*x1^2-0.167*x1*x3-1.16*x1-0.193*x3;
% -1.0*x3^2-0.397e-1*x3;
% -0.143*x1*x7-1.0*x7^2-0.388*x8*x7-0.163*x1-1.14*x7-0.443*x8;
% -1.0*x8^2-0.198*x8;
% -0.146*x9*x3-1.0*x9^2-0.136e-1*x9*x12-0.107*x3-0.730*x9-0.994e-2*x12;]

% [-1.0*x2^2-0.969*x4*x2-0.264*x5*x2-0.113*x10*x2-1.26*x2-1.22*x4-0.333*x5-0.142*x10;
% -0.866*x6*x4-1.0*x6^2-0.271*x4-0.313*x6;
%  -0.419*x6*x10-1.0*x10^2-0.822*x6-1.96*x10;
%  -0.548*x15*x3-0.796*x15*x6-0.230*x15*x10-1.0*x15^2-0.936*x3-1.36*x6-0.394*x10-1.71*x15;
%  -0.407*x16*x5-0.527*x15*x16-1.0*x16^2-0.275*x5-0.356*x15-0.676*x16];


% [-1.0*x4^2-0.926e-1*x4;
% -0.847*x4*x5-1.0*x5^2-1.77*x4-2.09*x5;
% -0.947e-1*x11*x10-1.0*x11^2-0.370e-2*x10-0.391e-1*x11;
% -0.801*x1*x12-1.0*x12^2-0.195*x1-0.244*x12;
% -0.211*x5*x13-1.0*x13^2-0.193*x5-0.915*x13;
% -0.571*x14*x9-0.709*x14*x12-1.0*x14^2-0.881*x9-1.09*x12-1.54*x14;
% ]



