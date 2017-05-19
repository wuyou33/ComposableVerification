function composable()
% parameters
direct_eigen_cal_on = false;

n=5;
blk_size=1;
num_blks=n/blk_size;

A=randn(n,n);



% density=.3;
% A=sprand(n,n,density)-10*eye(n);



if direct_eigen_cal_on
    direct_eigen_cal()
end

cvx_status='s';
while (~strcmp(cvx_status,'Solved'))
    cvx_begin sdp 
    cvx_solver Mosek
    variable P_i(blk_size,blk_size,num_blks) hermitian semidefinite
    tem_Cell=mat2cell(P_i,[blk_size],[blk_size],ones(1,num_blks));
    blkd_P=blkdiag(tem_Cell{:});
    subject to
    blkd_P >= 1e-6*eye(n)
    blkd_P*A+A'*blkd_P<= -1e-9*eye(n)
    cvx_end
end
lumped_P=full(blkd_P);
end



function Ricaati(A,n,blk_size)
scaling_choice = 'identity';
scaling_choice = 'weight';

if (strcmp(scaling_choice,'identity'))
    M=cell{ones};
else
    M=0;
end

%     A'XE + E'XA - (E'XB + S)R  (B'XE + S') + Q = 0
% When omitted, R, S and E are set to the default values R=I, S=0,
%     and E=I.
tic
[X,L,G,ricatti_flag]=care(A,A12,Q,-Q);
ricatti_time=toc
if ricatti_flag==-1 || ricatti_flag== -2
    ricatti_fails=true;
else
    ricatti_time;
end

end


function testing_arrows(A)
n=size(A,1);
for i =1:n
    difference=A(i,i);
    for j=1:n
        difference=difference-A(i,j)/A(j,j)*A(j,i);
        disp(A(i,i)-A(i,j)/A(j,j)*A(j,i));
    end
    disp('row difference');
    disp(difference)
end
end

function direct_eigen_cal(A)
    tic;
    find(real(eig(A))<0);
    toc
end

