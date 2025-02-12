function midx = Multi_Index(Nvar,order)
Nbasis = factorial(Nvar+order)/factorial(Nvar)/factorial(order);
midx = zeros(Nbasis,Nvar);
a = 1;
for i = 0:order
    if i == 0
        midx(1,:) = zeros(1,Nvar);
    else
        num = (factorial(Nvar + i-1)/factorial(Nvar)/factorial(i-1))*(Nvar/i);
        midx(a+1:a+num,:) = multi_idx(i,Nvar);
        a = a + num;
    end
end