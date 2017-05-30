
     % Compute the Lovasz number of the Petersen graph
     n = 10; % number of vertices
     % edges is a |E|x2 matrix
     edges = [[1,2];[2,3];[3,4];[4,5];[5,1];
              [1,6];[2,7];[3,8];[4,9];[5,10];
              [6,8];[6,9];[7,9];[7,10];[8,10]];
A = zeros(n,n);
     A(sub2ind(size(A),edges(:,1),edges(:,2))) = 1;
     A = A+A';
     cvx_begin sdp
    cvx_solver Mosek
    variable X(n,n) hermitian semidefinite
minimize(trace(A*(ones(n,n)-X))/4)
    subject to
    X >= 1e-9*eye(n)
    diag(X) == ones(n,1)
    cvx_end
    
   
%      X = sdpvar(n);
% %      Constraints = [X >= 0; diag(X) == ones(n,1)];
%      objective = -trace(A*(ones(n,n)-X))/4;
%      maxcutsol = solvesdp(constraints,objective);
%      double(objective)
