function [new_basis, E, E1, E2, E3, E_L2, init_TDPC] = updategPC_GramSchmidt(t,rv,Xe,xi_sample,o,Lb,Ld,TDbasis)


Nsample = length(TDbasis);
P = size(TDbasis,2);
x   =  Xe(1:P);
y   =  Xe(P+1:2*P);
th  =  Xe(2*P+1:3*P);
dx   =  Xe(3*P+1:4*P);
dy   =  Xe(4*P+1:5*P);
dth  =  Xe(5*P+1:6*P);

mean = sum(rv)/Nsample;
rv = rv - mean;

rvX = [TDbasis*x', TDbasis*y', TDbasis*th',TDbasis*dx', TDbasis*dy', TDbasis*dth'];

%% property
Nvar = size(rv,2);
order = o;
Nbasis = factorial(Nvar + order)/factorial(Nvar)/factorial(order);


%% multi index
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

%% multivariate monomials
Xbeta = zeros(Nsample, Nbasis);
for i = 1:Nbasis
    Pe = ones(Nsample,1);
    for k = 1:Nvar
        Pee = rv(:,k).^midx(i,k);
        Pe = Pe .* Pee;
    end
    Xbeta(:,i) = Pe;
end

%% Gram-Schmidt
new_basis = zeros(Nsample,Nbasis);
new_basis(:,1) = 1;
U_sample = zeros(Nbasis,Nbasis);
L_sample = zeros(Nbasis,Nbasis);
for i = 2:Nbasis
    new_basis(:,i) = Xbeta(:,i);
    for j = 1:i-1
        U = sum(Xbeta(:,i) .* new_basis(:,j))/Nsample;
        L = sum(new_basis(:,j).^2 )/Nsample;
        U_sample(i,j) = U;
        L_sample(i,j) = L;
        new_basis(:,i) = new_basis(:,i) - new_basis(:,j)*(U/L);
    end
end
% Normalize
for i = 1:Nbasis
    new_basis(:,i) = new_basis(:,i) / sqrt(sum(new_basis(:,i) .* new_basis(:,i))/Nsample);
end
%% modified Gram-Schmidt 
% new_basis = zeros(Nsample,Nbasis);
% 
% new_basis(:,1) = Xbeta(:,1) / sqrt(sum(Xbeta(:,1))/Nsample);
% for i = 2:Nbasis
%     new_basis(:,i) = Xbeta(:,i);
%     for j = 1:i-1
%         new_basis(:,i) = new_basis(:,i) - ((new_basis(:,j)'*new_basis(:,i))/Nsample) * new_basis(:,j);
%     end
%     new_basis(:,i) = new_basis(:,i) / sqrt(sum(new_basis(:,i).^2)/Nsample);
% end

%% Whitening transformation
G = zeros(Nbasis,Nbasis);
for i = 1:Nbasis
    for j = i:Nbasis
        Mij = sum(new_basis(:,i) .* new_basis(:,j))/Nsample;
        G(i,j) = Mij;
        G(j,i) = Mij;
    end
end
W = chol(G);
invW = W\eye(length(W));
psi = zeros(Nsample,Nbasis);
for i = 1:Nbasis
    for j = 1:i
        psi(:,i) = psi(:,i) + invW(j,i) * new_basis(:,j);
    end
end
new_basis = psi;

%% Galerkin projection
E = zeros(Nbasis,Nbasis);
for i = 1:Nbasis
    for j = i:Nbasis
        E(i,j) = sum(new_basis(:,i) .* new_basis(:,j))/Nsample;
        E(j,i) = E(i,j);
    end
end

E_L2 = zeros(Nbasis,Nbasis);
for i = 1:Nbasis
    for j = 1:Nbasis
        Mij = sum(new_basis(:,i) .* new_basis(:,j) .* (Lb + Ld * xi_sample).^2)/Nsample;
        E_L2(i,j) = Mij;
        E_L2(j,i) = Mij;
    end
end

E3 = zeros(Nbasis,1);
for i = 1:Nbasis
    Pe = sum(new_basis(:,i))/Nsample;
    E3(i) = Pe;
end

E1 = zeros(Nbasis,Nbasis^2);
for i = 1:Nbasis
    E1e = zeros(Nbasis,Nbasis);
    for j = 1:Nbasis
        for l = j:Nbasis
            Mij = sum(new_basis(:,i) .* new_basis(:,j) .* new_basis(:,l))/Nsample;
            E1e(j,l) = Mij;
            E1e(l,j) = Mij;
         end
    end
    E1(:,Nbasis*(i-1)+1 : Nbasis*i) = E1e;
end

E2 = zeros(Nbasis,Nbasis^3);
for i = 1:Nbasis
    for ii = 1:Nbasis
        ME2e = zeros(Nbasis,Nbasis);
        for j = 1:Nbasis
            for l = j:Nbasis
                Mij = sum(new_basis(:,i) .* new_basis(:,ii) .* new_basis(:,j) .* new_basis(:,l))/Nsample;
                ME2e(j,l) = Mij;
                ME2e(l,j) = Mij;
             end
        end
        n1 = Nbasis^2*(i-1) + Nbasis*(ii-1) + 1;
        n2 = Nbasis^2*(i-1) + Nbasis*(ii);
        E2(:,n1:n2) = ME2e;
    end
end

%% update Initial Condition

init_TDPC = zeros(Nbasis,6);
for j = 1:Nbasis
    
    init_TDPC(j,1) = sum( (Lb + Ld * xi_sample) .* cos(rvX(:,3)) .* new_basis(:,j))/Nsample;
    init_TDPC(j,2) = sum( (Lb + Ld * xi_sample) .* sin(rvX(:,3)) .* new_basis(:,j))/Nsample;
    init_TDPC(j,3) = sum(rvX(:,3) .* new_basis(:,j))/Nsample;
    
    init_TDPC(j,4) = sum( -(Lb + Ld * xi_sample) .* sin(rvX(:,3)) .* rvX(:,6) .* new_basis(:,j))/Nsample;
    init_TDPC(j,5) = sum(  (Lb + Ld * xi_sample) .* cos(rvX(:,3)) .* rvX(:,6) .* new_basis(:,j))/Nsample;
    init_TDPC(j,6) = sum(rvX(:,6) .* new_basis(:,j))/Nsample;
    
end

init_TDPC = reshape(init_TDPC,Nbasis*6,1);

end