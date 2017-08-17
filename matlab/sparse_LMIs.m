function sparse_LMIs(A,n,blk_size)
load('all_A.mat');
n=size(all_A,1);
blk_size=size(all_A,2);
num_samples=size(all_A,3);
num_blks=n/blk_size;
samplecount=0;

cvx_status='s';
cvx_begin sdp quiet
cvx_solver Mosek



while (sample_counter<=num_samples)
    
    variable P_i(blk_size,blk_size,num_blks) hermitian semidefinite;
    variable M_ij(blk_size,blk_size,num_blks*(num_blks-1)) hermitian semidefinite;
    
    
    tem_Cell_M=mat2cell(M_ij,[blk_size],[blk_size],ones(1,num_blks*(num_blks-1)));
    % blkd_M=blkdiag(tem_Cell_M{:});
    
    array_of_D=zeros(n,n,num_blks);
    
    for i =1:num_blks
        for j = 1:n
            if i==j
            	array_of_D(j*blk_size:(j+1)*blk_size,:,i)=P_i;
            	array_of_D(:,j*blk_size:(j+1)*blk_size,i)=P_i;
            	array_of_D(j*blk_size:(j+1)*blk_size,j*blk_size:
            		(j+1)*blk_size,i)=P_i*A_ii+A_ii'*P_i;
            else
            	array_of_D(j*blk_size:(j+1)*blk_size,j*blk_size:
            		(j+1)*blk_size,i)=-tem_Cell_M(i*j);
            end            
        end
    end
end

tem_Cell_D=mat2cell(D,[n],[n],ones(1,num_blks));
blkd_M=blkdiag(tem_Cell_D{:});

subject to
% >= 1e-7*eye(n)
blkd_P*A+A'*blkd_P<= -1e-9*eye(n)
cvx_end

end
