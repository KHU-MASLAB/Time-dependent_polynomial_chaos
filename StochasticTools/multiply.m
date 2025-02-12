function EA = multiply(E,A,dim)

P = size(E,1);
D1 = size(A,1)/(P^dim);
D2 = size(A,2);
EA = zeros(D1*P,D2*P);
for i = 1:P
    n1 = (i-1)*(P^dim) + 1;
    n2 = i*(P^dim);
    EE = E(:,n1:n2);
    for j = 1:D1-1
        EE = blkdiag(EE,E(:,n1:n2));
    end   
    for k = 1:D2
        EA(:,i+(k-1)*P) = EE*A(:,k);
    end    
end