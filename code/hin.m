not_stable = true;

while not_stable

A11=randn(3);
A22=randn(1);
if real(eig(A11))<0
    if real(eig(A22))<0
        not_stable = false;
    end;
end;
end;


% P1=lyap(A11', eye(3));
% P2=lyap(A22', eye(1));

disp('the rho part starts')

% new loop
cvx_begin sdp
variable rho(1);

A12=[0;0;rho];
A13=[0;0;0.01];
A21=[0,0,rho];

% A=[A11,A12;A21,A22];

variable P1(3,3) symmetric;
% P1 = semidefinite(3);

variable P2(1,1) symmetric;
% P2 = semidefinite(1);

% P = [P1,zeros(3,1);zeros(1,3),P2];

% Dot=P*A+A'*P;
minimize(-rho)
subject to
[P1*A11+A11'*P1+eye(3),P1*A12;A12'*P1,eye(3)]<=zeros(3,3);
% P1*A11+A11'*P1+P1*A13*A13'*P1+eye(2)<=zeros(3,3);
% 	Dot <= -(1e-3)*eye(4);
cvx_end
disp(norm(inv(A11)*A12));
disp(norm(inv(A22)*A21));



% % testing if can do better by changing the 
