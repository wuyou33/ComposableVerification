function [ricatti_succ_count, ricatti_time] = Ricaati_eqns()
n=10;
blk_size=2;
num_blks=n/blk_size;
num_samples=50;

Asize=num2str(n);
ABlkSize=num2str(blk_size);
load(strcat('A_','size_',Asize,'BlkSize_',ABlkSize,'num_samples_',num_samples,'.mat'))

scaling_choice = 'identity';
% scaling_choice = 'sigma1';
% When omitted, R, S and E are set to the default values S=0,
% and E=I. So setting  R=-I, solving X for A'X + XA - (XB)(B'X) + Q = 0
% or equivalently solving for dual PA'+AP+BB'+PQP=0, setting Q=eye then.
Q=eye(blk_size);
R=-eye(blk_size);
S=zeros(blk_size);
E=Q;

ricatti_time=0;
ricatti_succ_count=0;
for k=1:num_samples
    A=all_A{k};
    M=scaling_matrices(n,blk_size,num_blks,A,scaling_choice);
     for i=1:num_blks
        Aii=A((i-1)*blk_size+1:(i)*blk_size,(i-1)*blk_size+1:(i)*blk_size);
        A((i-1)*blk_size+1:(i)*blk_size,(i-1)*blk_size+1:(i)*blk_size)=zeros(blk_size,blk_size);
        tic;
        [X,L,G,ricatti_flag]=care(Aii,A([(i-1)*blk_size+1:(i)*blk_size],:)*[M{i,:}]',Q,R,S,E);
        ricatti_time_incremental=toc;
        if ricatti_flag == -1 || ricatti_flag == -2
            break
        end
        ricatti_time=ricatti_time+ricatti_time_incremental;
    end
    if ricatti_flag ~= -1 & ricatti_flag ~= -2
        ricatti_succ_count = ricatti_succ_count+1;
    end
end
ricatti_time=ricatti_time/num_samples;
end

function M=scaling_matrices(n,blk_size,num_blks,A,scaling_choice)
    if (strcmp(scaling_choice,'identity'))
        M=repmat(ones(blk_size,blk_size),n/blk_size,n/blk_size);
        M=mat2cell(M,blk_size*ones(1,num_blks),blk_size*ones(1,num_blks));
    else
        A=mat2cell(A,blk_size*ones(1,num_blks),blk_size*ones(1,num_blks));
        % use the largest singular value as a huristic guess
        r=zeros(num_blks,num_blks);
        M={};
        for i=1:num_blks
            for j=1:num_blks
                r(i,j)=svds(inv(A{i,i})*A{i,j},1);
            end
        end

        for i=1:num_blks
            for j=1:num_blks
                if r(i,j)*r(j,i)>1 and r(i,j)>1
                    r(i,j)=1/r(j,i);
                end
                M{i,j}=sqrt(r(i,j))*eye(blk_size);
            end
        end
    end                
    for i=1:num_blks
        M{i,i}=zeros(blk_size,blk_size);
    end
end


