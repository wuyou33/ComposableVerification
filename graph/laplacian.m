L = laplacian(A);
>> [V D] = eigs(L, 2, 'SA');
>> D(2,2)

ans =

    0.0738
    
>> plot(sort(V(:,2)), '.-');
>> [ignore p] = sort(V(:,2));
>> spy(A(p,p));