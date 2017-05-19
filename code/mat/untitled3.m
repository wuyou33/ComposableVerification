not_stable = true;

while not_stable

A11=randn(3);
A22=randn(3);
if real(eig(A11))<0
    if norm(A11)>1
        if real(eig(A22))<0
            if norm(A22)>1
                not_stable = false;
            end;
        end;
    end;
end;
end;

A12=zeros(3,3);
A21=zeros(3,3);


A12(3)=1e-1;
A21(3)=1e-1;
A13=zeros(3,3);
A13(2)=2e-1;

% P1=lyap(A11', eye(3));
% P2=lyap(A22', eye(1));

% disp('the rho part starts')

% new loop
cvx_begin quiet sdp
cvx_solver mosek
variable gamma12;
variable gamma13;
variable gamma21;
% A12=[0;0;rho];
% A13=[0;0;0.01];
% A21=[0,0,rho];

% A=[A11,A12;A21,A22];

% variable P1(3,3) symmetric;
P1 = semidefinite(3);

% variable P2(3,3) symmetric;
P2 = semidefinite(3);
P=blkdiag(P1,P2)
% P = [P1,zeros(3,1);zeros(1,3),P2];

% Dot=P*A+A'*P;
% minimize(gamma12+gamma13)
subject to
P>=eye(6)
% [A11*Q1+Q1*A11'+eye(3),P1*A12;A12'*P1,eye(3)]<=zeros(3,3);
% P1*A11+A11'*P1+P1*A13*A13'*P1+eye(2)<=zeros(3,3);

% [P1*A11+A11'*P1+eye(3),P1*A12;(P1*A12)',-gamma12*eye(3)]<0;
% [P1*A11+A11'*P1+eye(3),P1*A13;(P1*A13)',-gamma13*eye(3)]<0;

% [P2*A22+A22'*P2+eye(3),P2*A21;(P2*A21)',-gamma21.*eye(3)]<=0;

cvx_end

% Q1=inv(P1);
% Q2=inv(P2);
% ricatti1=P1*A11+A11'*P1+eye(3)+P1*A12*(P1*A12)'/gamma1
% gamma1;
% gamma2;

% norm(inv(A11*Q1+Q1*A11')*(A12*Q2+Q1*A21'))


% [A11*Q1+Q1*A11',A12*Q2+Q1*A21';(A12*Q2+Q1*A21')',A22*Q2+Q2*A22' ]
% disp(norm(inv(A11)*A12));
% disp(norm(inv(A22)*A21));



% % testing if can do better by changing the 
