function [ricatti_succ_count, avg_Riccati_time] = Ricaati_eqns(sample_A)


n=10;
blk_size=2;
num_blks=n/blk_size;
num_samples=10;

% scaling_choice = 'identity';
scaling_choice = 'sigma1';

% if from_file
%     Asize=num2str(n);
%     ABlkSize=num2str(blk_size);
%     load(strcat('A_','size_',Asize,'BlkSize_',ABlkSize,'num_samples_',num_samples,'.mat'),'sample_A');
% end

% R,S,E will be used later in the call of the built-in Riccati solver care()
% Check the docmentation of care() for term matching
% With the following setup: R=-I,S=0,E=I, we are solving X for 
% A'X + XA + (XB)(B'X) + Q = 0
R=-eye(blk_size);
S=zeros(blk_size);
E=eye(blk_size);
% where in our problem setting, A corresponds to A_{ii}',B corresponds to the square root of the sum of {inv{M_{ij}} over j, and Q corresponds to the sum of all the A_{ij}M_{ij}A_{ij}' terms

total_time=0;
ricatti_succ_count=0;
for k=1:num_samples
    A=sample_A{k};
    A=mat2cell(A,blk_size*ones(1,num_blks),blk_size*ones(1,num_blks));
    M=scaling_matrices(n,blk_size,A,scaling_choice);
    for i=1:num_blks
        this_sample_time=0;
        Q=zeros(blk_size,blk_size);
        B=zeros(blk_size,blk_size);
        % Aii=A((i-1)*blk_size+1:(i)*blk_size,(i-1)*blk_size+1:(i)*blk_size);
        for j=1:num_blks
            if i~=j
                Q=Q+A{i,j}*M{i,j}*A{i,j}';
                B=B+inv(M{i,j}); 
            end
        end
        B=chol(B)';
        tic;
        [X,L,G,ricatti_flag]=care(A{i,i}',B,Q,R,S,E);
        this_block_time=toc;
        if ricatti_flag == -1 || ricatti_flag == -2
            % -1 if the Hamiltonian matrix has jw-axis eigenvalues
            % -2 if there is no finite stabilizing solution X
            % meaning the solve is not successful
            ricatti_succ_flag=0;
            this_sample_time=0;
            break
        else
            % solve is successful, update time and count
            ricatti_succ_flag=1;
            this_sample_time=this_sample_time+this_block_time;
        end
    end
    ricatti_succ_count = ricatti_succ_count+ricatti_succ_flag;
    total_time=total_time+this_sample_time;
end

avg_Riccati_time=total_time/ricatti_succ_count;

end

function M=scaling_matrices(n,blk_size,A,scaling_choice)
num_blks=n/blk_size;
M=repmat(eye(blk_size,blk_size),n/blk_size,n/blk_size);
M=mat2cell(M,blk_size*ones(1,num_blks),blk_size*ones(1,num_blks));

if (strcmp(scaling_choice,'sigma1'))
    % A=mat2cell(A,blk_size*ones(1,num_blks),blk_size*ones(1,num_blks));
    % use the largest singular value as a huristic guess
    r=zeros(num_blks,num_blks);
    for i=1:num_blks
        for j=1:num_blks
            r(i,j)=svds(inv(A{i,i})*A{i,j},1);
            % disp(r(i,j))
        end
    end
    
    for i=1:num_blks
        for j=1:num_blks
            if r(i,j)*r(j,i)>1
                if r(i,j)>1
                    r(i,j)=1/r(j,i);
                    M{i,j}=sqrt(r(i,j))*M{i,j};
                end
            end
            % M{i,j}=sqrt(r(i,j))*M{i,j};
            % disp(M{i,j})
        end
    end
end

end


