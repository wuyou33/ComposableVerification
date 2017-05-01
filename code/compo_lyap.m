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


P1=lyap(A11', eye(3));
P2=lyap(A22', eye(1));


disp('the rho part starts')

% new loop
cvx_begin sdp
variable rho(1);

A12=[0;0;rho];
A21=[0,0,rho];

A=[A11,A12;A21,A22];

% variable P1(2,2) symmetric;
% P1 = semidefinite(2);

% variable P2(1,1) symmetric;
% P2 = semidefinite(1);

P = [P1,zeros(3,1);zeros(1,3),P2];

Dot=P*A+A'*P;
minimize(-rho)
subject to
	Dot <= -(1e-3)*eye(4);
cvx_end
disp(norm(inv(A11)*A12));
disp(norm(inv(A22)*A21));



% % testing if can do better by changing the 

function testing_scaling(A11,A12,A21,A22)

end 