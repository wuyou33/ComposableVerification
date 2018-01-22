blk_size=2;
num_blks=5;
Mij=randn(blk_size,blk_size,num_blks*(num_blks-1));
M_ij=Mij+permute(Mij,[2 1 3]);

disp(M_ij);
M_ij=permute(M_ij,[3,4,1,2]);
    M_ij=reshape(M_ij,num_blks-1,num_blks,blk_size,blk_size);
    M_ij(:,num_blks+1,:,:)=zeros(size(M_ij(:,1,:,:)));
    M_ij=permute(M_ij,[2,1,3,4]);
    M_ij=reshape(M_ij,[(num_blks-1)*(num_blks+1),blk_size,blk_size,]);
    M=cat(1,zeros(size(M_ij(1,:,:))),M_ij);
    M=reshape(M,num_blks,num_blks,blk_size,blk_size);
    M=permute(M,[3,4,1,2]);
