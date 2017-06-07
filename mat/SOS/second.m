degV=4;

x=msspoly('x',12);
y=msspoly('y',2);

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
y1=y(1);
y2=y(2);

X1 = [x1, x2, x7, x9, x11, x12]';

X2 = [x1,x3, x4, x5, x6, x8, x10]';

dx1dt = -1.57*x1 - 0.574*x2 + 0.861*x7  + 1.10*x11 + 1.14*x12+ 1.51*x9 ;- 0.583*x3 + 0.942*x4 + 1.35*x6  +x2*x9*x6
dx2dt = -0.760*x2 - 1.06*x7  + 0.681*x11 ;+1.54*x3 - 1.63*x4 - 0.819*x5 - 0.925*x6 + x7*x10*x11- 0.535*x10;
dx7dt =-0.845*x1- 1.67*x7 + 1.72*x9   - 0.600*x12 + x1*x7*x11- 1.25*x11;+ 0.839*x10
dx9dt = -1.25*x1 + 1.04*x2  - 1.20*x9 - 3.15*x11 - 0.998*x12 ;+ 0.990*x3 - 1.44*x4 - 1.44*x5 -1.29*x6- 1.16*x10 + x1*x6*x10
dx11dt = + 2.00*x7  + 2.46*x9  - 0.833*x11 + 1.29*x12 ;-2.29*x3 + 1.41*x4 + 1.43*x5 + 0.512*x6 - 0.801*x8- 0.659*x10+ x6^3
dx12dt = -0.850*x1 + 0.515*x7 + 1.20*x9 - 1.25*x11 + x1*x9*x11;- 0.919*x3 
dX1dt=[dx1dt;dx2dt;dx7dt;dx9dt;dx11dt;dx12dt];

dx1dt = -1.57*x1 - 0.583*x3 + 0.942*x4 + 1.35*x6; - 0.574*x2+ 0.861*x7  + 1.51*x9 + 1.10*x11 + 1.14*x12 +x6*x2*x9
dx3dt =  1.77*x1 - 0.983*x2 - 0.766*x3 - 0.930*x6 - 1.12*x9 - 0.998*x10 + 1.84*x11 + 0.772*x12 +x3^3;
dx4dt = -1.37*x1 + 1.31*x2 + 0.790*x3 - 0.572*x4 - 1.07*x5 - 0.783*x6 - 0.938*x8 + 2.03*x9 - 0.857*x11 + x4*x8*x11; 
dx5dt = -1.04*x1 + 0.755*x2 + 1.12*x4 - 0.538*x5 - 0.563*x8 + 1.40*x9 - 1.32*x11 + x1*x2*x5;
dx6dt =  1.26*x2 + 1.11*x4 - 0.885*x6 - 0.935*x8 + 1.05*x9 - 1.30*x11 + x2*x4*x11;
dx8dt =  0.583*x1 + 1.20*x4 + 0.816*x5 + 0.599*x6 - 0.735*x9 + x1^3;
dx10dt =  1.43*x2 + 0.835*x3 + 0.545*x7 - 0.890*x8 - 0.973*x10 + x2*x3*x8;
dX2dt=[dx1dt;dx3dt;dx4dt;dx5dt; dx6dt;dx8dt;dx10dt];



% construct V
V1Xmonom = monomials(X1,0:degV);
V2Xmonom = monomials(X2,0:degV);

% partition 1
prog = spotsosprog;
prog = prog.withIndeterminate([X1;y1]);

%V condition
[prog,V1] = prog.newFreePoly(V1Xmonom);
prog = prog.withSOS(V1);
PV1PX1 = diff(V1,X1);
V1dot = diff(V1,X1)*dX1dt;
% V1dotcomp=-V1dot+2*X1'*PV1PX1'*(0.136e-1*x9)*y(1)+y(1)^2-(-0.936*x3)^2-(-0.195*x1)^2-(x1)^2-(-0.881*x9)^2-(x9)^2;
prog = prog.withSOS(-V1dot);
% prog = prog.withSOS(-V1dotcomp);

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

disp(V1);
% disp(V2);





% V1(X1) = 1.59*x8*x4 - 1.52*x5*x4 + 1.35*x3^2 + 2.03*x10^2 + 2.33*x8^2 + 2.42*x5^2 + 2.17*x6^2 + 3.77*x4^2 +  2*x1^2 -0.899*x5*x8 
% + 0.993*x4*x3 - 0.181*x4*x6 - 1.08*x4*x10 + 1.04*x8*x1 + 0.182*x8*x3 + 0.627*x8*x6 - 0.250*x8*x10 - 0.641*x5*x3 + 0.793*x5*x6 
% + 0.621*x5*x10 - 0.321*x10*x3 - 0.432*x10*x6 + 0.652*x1*x3 + 0.989*x1*x4 - 0.755*x1*x10 - 0.286*x1*x5 - 0.334*x1*x6-.524*x6*x3
% V2(X2) = 1.86*x9^2 + 0.355*x9*x2 + 2.09*x2^2 + 0.0663*x1*x9 -0.226*x1*x2 + 1.49*x1^2 - 0.04751*x9*x11 + 0.101*x2*x11  
% + 0.383*x1*x11 + 2.16*x11^2 + 0.602*x9*x7 - 0.363*x2*x7 + 0.0977*x1*x7 + 0.216*x11*x7 + 1.50*x7^2 + 0.143*x9*x12
% - 0.0914*x12*x2 - 0.261*x1*x12 - 0.388*x12*x11 + 0.0748*x12*x7 + 2.23*x12^2

% V2_sparse(X2) = 0.101*x9^2 + 0.1*x2^2 - 0.421e-7*x1*x9 + 0.1*x1^2 - 0.445e-3*x11*x9 + 0.1*x11^2 + 0.12*x7^2 + 0.277e-2*x12*x9 
% - 0.121e-2*x12*x11 + 0.372e-7*x12*x7 + 0.102*x12^2