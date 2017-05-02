function Ginf=Hinfnorm(A,B,C,D)
dim=size(B);
n=dim(1);
m=dim(2);
cvx_begin sdp
    variable P(n,n) symmetric
    variable g;
    minimize g;
    
    subject to
        P>=0;
        [A'*P+P*A+C'*C  , P*B+C'*D;
        (P*B+C'*D)'    D'*D-g*eye(m)]<=0;
Ginf=sqrt(g);