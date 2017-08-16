function [ricatti_succ_count, ricatti_time] = Ricaati_eqns(A,n,blk_size,scaling_choice)
scaling_choice = 'identity';
scaling_choice = 'sigma';
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
% copy_total=A_total_cell;

for i=1:num_blks
    A([(i-1)*blk_size+1:(i)*blk_size],[(i-1)*blk_size+1:(i)*blk_size]) =zeros(blk_size,blk_size);
    % copy_total=mat2cell(copy_total,blk_size*ones(1,num_blks),blk_size*(num_blks-1))
end

for i=1:num_blks
    tic
    [X,L,G,ricatti_flag]=care(A_total_cell{i,i},A([(i-1)*blk_size+1:(i)*blk_size],:),Q,-R);
    ricatti_time_incremental=toc
    ricatti_time=ricatti_time+ricatti_time_incremental;
end
if ricatti_flag == -1 || ricatti_flag == -2
    ricatti_succ_count = 0;
else
    ricatti_succ_count = 1;
end
end
