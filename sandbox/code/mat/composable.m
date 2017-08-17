function [ricatti_count, ricatti_time, lumped_P]=composable()
% parameters
direct_eigen_cal_on = false;
A_option = 'sparse_A';
% A_option = 'dense_A';

% the LTI example from James Anderson's paper but seems wrong as 
% A=[-2 4 0 0;-2 -2 -4 -6;-2 2 -2 0;-2 2 4 -2];

A=[-2 0 4 0;-2 -2 2 0;-2 -4 -2 -6;-2 2 4 -2];

n=4;
blk_size=2;
num_samples=1000;

num_blks=n/blk_size;
% if strcmp(A_option, 'sparse_A')
%     density=.2;
%     A=sprand(n,n,density)-10*eye(n);
% elseif strcmp(A_option, 'dense_A')
%     A=randn(n,n);
% else ;
% end

A=[-2 4 0 0;-2 -2 -4 -6;-2 2 -2 0;-2 2 4 -2];

if direct_eigen_cal_on
    direct_eigen_cal();
end

cvx_status='s';
sample_counter=0;
ricatti_count=0;

while (sample_counter<=num_samples)
    cvx_begin sdp
    cvx_solver Mosek
    variable P_i(blk_size,blk_size,num_blks) hermitian semidefinite
    tem_Cell=mat2cell(P_i,[blk_size],[blk_size],ones(1,num_blks));
    blkd_P=blkdiag(tem_Cell{:});
    subject to
    blkd_P >= 1e-7*eye(n)
    blkd_P*A+A'*blkd_P<= -1e-9*eye(n)
    cvx_end
    if (strcmp(cvx_status,'Solved'))% the lumped version is successful
        sample_counter=sample_counter+1;
        [ricatti_succ_count, ricatti_time] = Ricaati(A,n,blk_size,num_blks)
%         ricatti_count=ricatti_count+ricatti_succ_count;
    end
end
lumped_P=full(blkd_P);
end

function sparse_LMIs(A,n,blk_size)
cvx_status='s';
cvx_begin sdp
cvx_solver Mosek
variable P_i(blk_size,blk_size,num_blks) hermitian semidefinite;
tem_Cell_P=mat2cell(P_i,[blk_size],[blk_size],ones(1,num_blks));
blkd_P=blkdiag(tem_Cell_P{:});
P_zero=mat2cell(zeros(n,n,n),n,n,ones(1,n));
scaling_zero=0;
% P_zero(:,:,) = tem_Cell_P(i);
variable M_i(blk_size,blk_size,num_blks*(num_blks-1)) hermitian semidefinite;
tem_Cell_M=mat2cell(M_i,[blk_size],[blk_size],ones(1,num_blks*(num_blks-1)));
blkd_M=blkdiag(tem_Cell_M{:});
subject to
 % >= 1e-7*eye(n)
blkd_P*A+A'*blkd_P<= -1e-9*eye(n)
cvx_end
end


function [ricatti_succ, ricatti_time] = Ricaati(A,n,blk_size,num_blks)
scaling_choice = 'identity';
scaling_choice = 'weight';

if (strcmp(scaling_choice,'identity'))
    M=cell{ones};
else
    % use the largest singular value as a huristic guess
    M=svds(A,1);
end


Q=num_blks*eye(blk_size);

R=num_blks*eye(n);
% A'XE + E'XA - (E'XB + S)R  (B'XE + S') + Q = 0
% When omitted, R, S and E are set to the default values R=I, S=0,
% and E=I.
ricatti_time=0;
A_total_cell=mat2cell(A,blk_size*ones(1,num_blks),blk_size*ones(1,num_blks));
copy_total=A_total_cell;

for i=1:num_blks
    A([(i-1)*blk_size+1:(i)*blk_size],[(i-1)*blk_size+1:(i)*blk_size]) =zeros(blk_size,blk_size);
    % copy_total=mat2cell(copy_total,blk_size*ones(1,num_blks),blk_size*(num_blks-1))
end

ricatti_succ_count=0;

for i=1:num_blks
    tic
    [X,L,G,ricatti_flag]=care(A_total_cell{i,i},A([(i-1)*blk_size+1:(i)*blk_size],:),Q,-R);
    disp(ricatti_flag);
    ricatti_time_incremental=toc
    ricatti_time=ricatti_time+ricatti_time_incremental;
    if ricatti_flag~=-1 & ricatti_flag~= -2
        ricatti_succ_count=ricatti_succ_count+1;
    end
end

if ricatti_succ_count==num_blks
    disp('true')
    ricatti_succ=true;
else
    ricatti_succ=false;
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

