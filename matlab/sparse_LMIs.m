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
