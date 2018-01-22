function sparse_time=sparse_LMIs(sample_A)

% load('all_A.mat');
% n=size(sample_A,1);
% blk_size=size(all_A,2);
% num_samples=size(all_A,3);
% num_blks=n/blk_size;
n=10;
blk_size=2;
num_blks=n/blk_size;
num_samples=10;
total_time=0;




for k=1:num_samples
    A=sample_A{k};
    A=mat2cell(A,blk_size*ones(1,num_blks),blk_size*num_blks);
    cvx_status='s';
    cvx_begin sdp quiet
    cvx_solver Mosek
    % for each sample, create the blk-diag P matrix
    variable P_i(blk_size,blk_size,num_blks) hermitian semidefinite;
    variable M_ij(blk_size,blk_size,num_blks*(num_blks-1)) hermitian semidefinite;

    P_Cell=mat2cell(P_i,[blk_size],[blk_size],ones(1,num_blks));
    D_Cell=mat2cell(zeros(n,n,num_blks),[blk_size*ones(1,num_blks)],[blk_size*ones(1,num_blks)],[ones(1,num_blks)]);
    % blkd_P=blkdiag(P_Cell{:});
    
    % what's this for loop for??
    % for i = 1:num_blks
    %     wew = zeros(n,n);
    %     wew(i*blk_size:(i+1)*blk_size,i*blk_size:
    %                 (i+1)*blk_size)=tem_Cell{:,:,i}
    % end
    % still for each sample, create the auxilary M_ij matices
    M=permute(M_ij,[3,4,1,2]);
    M=reshape(M,num_blks-1,num_blks,blk_size,blk_size);
    M(:,num_blks+1,:,:)=zeros(size(M(:,1,:,:)));
    M=permute(M,[2,1,3,4]);
    M=reshape(M,[(num_blks-1)*(num_blks+1),blk_size,blk_size,]);
    M=cat(1,zeros(size(M(1,:,:))),M);
    M=reshape(M,num_blks,num_blks,blk_size,blk_size);
    M=permute(M,[3,4,1,2]);

    for i=1:num_blks
        % D_Cell{i,i,i}=D_Cell{i,i,i}+reshape(sum(reshape(M(:,:,i,:),num_blks,blk_size,blk_size),1),blk_size,blk_size);
        M_temp=mat2cell(M(:,:,i,:),blk_size,blk_size,1,ones(1,num_blks));
        D_Cell{i,i,i}=D_Cell{i,i,i}+sum(M(:,:,i,:),4);
        PA_row=P_Cell{i}*(A{i});
        D_Cell(i,:,i)=cellfun(@plus,D_Cell(i,:,i),mat2cell(PA_row,blk_size,blk_size*ones(1,num_blks)),'uni',0);
        D_Cell(:,i,i)=cellfun(@plus,D_Cell(:,i,i),mat2cell(PA_row',blk_size*ones(1,num_blks),blk_size),'uni',0);
        % D_Cell{:,i,i}=D_Cell{:,i,i}+P_Cell{i}*(A{i})';
        D_Cell(:,:,i)=cellfun(@minus,D_Cell(:,:,i),mat2cell(blkdiag(M_temp{:}),[blk_size*ones(1,num_blks)],[blk_size*ones(1,num_blks)]),'uni',0);
    end

    csize = size(D_Cell);
    % Treat 3+ dimension arrays

    % Construct the matrix by concatenating each dimension of the cell array into
    %   a temporary cell array, CT
    % The exterior loop iterates one time less than the number of dimensions,
    %   and the final dimension (dimension 1) concatenation occurs after the loops

    % Loop through the cell array dimensions in reverse order to perform the
    %   sequential concatenations
    for cdim=(length(csize)-1):-1:1
        % Pre-calculated outside the next loop for efficiency
        ct = cell([csize(1:cdim) 1]);
        cts = size(ct);
        ctsl = length(cts);
        mref = {};

        % Concatenate the dimension, (CDIM+1), at each element in the temporary cell
        %   array, CT
        for mind=1:prod(cts)
            [mref{1:ctsl}] = ind2sub(cts,mind);
            % Treat a size [N 1] array as size [N], since this is how the indices
            %   are found to calculate CT
            if ctsl==2 && cts(2)==1
                mref = {mref{1}};
            end
            % Perform the concatenation along the (CDIM+1) dimension
            ct{mref{:}} = cat(cdim+1,D_Cell{mref{:},:});
        end
        % Replace M with the new temporarily concatenated cell array, CT
        D_Cell = ct;
    end

    % Finally, concatenate the final rows of cells into a matrix
    m = cat(1,D_Cell{:});    
    % D=reshape(cell2mat(D_Cell),n,n*num_blks)
    % D=mat2cell(D,[n],n*ones(1,num_blks))
    % reshape(D_Cell,[n,n,num_blks])
    
    % tem_Cell_M=mat2cell(M_ij,[blk_size],[blk_size],ones(num_blks,(num_blks-1)));
    % add padding zero blocks
    
    % reshape the cell into a 1D cell
    %    for i = 1:num_blks
    %    	for j=1:num_blks
    %    		if i==j
    %     		MIJ(i,j)=zeros(blk_size,blk_size);
    %     	else
    %     		MIJ(i,j)=tem_Cell_M{i*j}
    %     	end
    %     end
    % end
    
    %    array_of_D=zeros(n,n,num_blks);
    
    %    for i =1:num_blks
    %    	A_ijM_ijAij=zeros(blk_size,blk_size)
    %        for j = 1:num_blks
    %        	A_ii=all_A(i*blk_size:(i+1)*blk_size,i*blk_size:
    %            		(i+1)*blk_size,i);
    %        	A_ij=all_A(i*blk_size:(i+1)*blk_size,j*blk_size:
    %            		(j+1)*blk_size,i);
    %            if i~=j
    %            	array_of_D(j*blk_size:(j+1)*blk_size,j*blk_size:
    %            		(j+1)*blk_size,i)=-tem_Cell_M(i*j);
    %            	A_ijM_ijAij=A_ijM_ijAij+A_ij*tem_Cell_M(i*j)*A_ij;
    %            else
    %            	array_of_D(j*blk_size:(j+1)*blk_size,:,i)=P_i;
    %            	array_of_D(:,j*blk_size:(j+1)*blk_size,i)=P_i;
    %            	array_of_D(j*blk_size:(j+1)*blk_size,j*blk_size:
    %            		(j+1)*blk_size,i)=-(P_i*A_ii+A_ii'*P_i+A_ijM_ijAij);
    %            end
    %        end
    %    end
    D=mat2cell(m,[n],[n],ones(1,num_blks));
    % blkd_M=blkdiag(tem_Cell_D{:});
    
    subject to
    blkdiag(D{:})<=-1e-7*eye(num_blks*n)
    blkdiag(P_Cell{:})>=1e-7*eye(n)
    % >= 1e-7*eye(n)
    % blkd_P*A+A'*blkd_P<= -1e-9*eye(n)
    tic
    cvx_end
    this_time=toc;
    total_time=total_time+this_time;
end
sparse_time=total_time/num_samples;
end