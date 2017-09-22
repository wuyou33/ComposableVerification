function [all_A, avg_lumped_time]=Random_A_Matrices(n,blk_size,num_samples)
% parameters
direct_eigen_cal_on = false;
A_option = 'sparse_A';
% A_option = 'dense_A';

num_blks=n/blk_size;
all_A=zeros(n,n,num_samples);
if strcmp(A_option, 'sparse_A')
    density=.2;
    A=sprand(n,n,density)-5*eye(n);
elseif strcmp(A_option, 'dense_A')
    A=randn(n,n);
end

cvx_status='s';
sample_counter=0;
lumped_total_time=0;
all_A={};
while (sample_counter<num_samples)
    tic
    cvx_begin sdp quiet
    cvx_solver Mosek
    variable P_i(blk_size,blk_size,num_blks) hermitian semidefinite
    tem_Cell=mat2cell(P_i,[blk_size],[blk_size],ones(1,num_blks));
    blkd_P=blkdiag(tem_Cell{:});
    subject to
    blkd_P >= 1e-7*eye(n)
    blkd_P*A+A'*blkd_P <= -1e-7*eye(n)
    cvx_end
    lump=toc;
    if (strcmp(cvx_status,'Solved'))% the lumped version is successful
        sample_counter=sample_counter+1;
        all_A(sample_counter)={A};
    end
    lumped_total_time=lumped_total_time+lump;
end
avg_lumped_time=lumped_total_time/num_samples;

Asize=num2str(n);
ABlkSize=num2str(blk_size);
save(strcat('A_','size_',Asize,'BlkSize_',ABlkSize,'num_samples_',num_samples,'.mat'),'all_A');
end


% testing direct eigenvalue based methods
function eigen_time=direct_eigen_cal(A)
tic;
find(real(eig(A))<0);
eigen_time=toc;
end

% function testing_arrows(A)
% n=size(A,1);
% for i =1:n
%     difference=A(i,i);
%     for j=1:n
%         difference=difference-A(i,j)/A(j,j)*A(j,i);
%         disp(A(i,i)-A(i,j)/A(j,j)*A(j,i));
%     end
%     disp('row difference');
%     disp(difference)
% end
% end


    