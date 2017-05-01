U1 = [100;1];
U2 = [1;100];
P1=diag(U1);
P2=diag(U2);

% % % looking for the touching hyperplane
cvx_begin 
variable x1(2);
% variable x2(2);

maximize ([10,0]*x1)

subject to 
[1,x1';x1,inv(P1)]>=1e-7*eye(3);
% [1,x2';x2,inv(P2)]>=1e-7*eye(3);
% P1*x1==P2*x2;

cvx_end