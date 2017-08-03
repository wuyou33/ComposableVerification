A =[-1.0000  0.2   0.7;
    0.2  -1.0000    0.7;
    0.9 0.09  -1.0000];
% % 
% 
% A =[-1.0000  0   ;
%     1.1  -1.0000];
  disp(real(eig(A)))  
% A=-eye(3)
cvx_begin
variable p1(1)
variable p2(1)
variable p3(1)

% variable p4(1)
% P=[p1, 0 p4;0 p2 0;p4 0 p3]
% minimize norm(p4)

P=[p1, 0 0;0 p2 0;0 0 p3]
% P=[p1,0;0,p2];
maximize p2
subject to 
p1>=1e-10;
p1<=1
p2>=1e-8;
p2<=1
p3>=0;
P*A+A'*P<=zeros(3,3);
cvx_end
disp(P)