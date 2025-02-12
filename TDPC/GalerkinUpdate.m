function Galerkin = GalerkinUpdate(basis,XiSample,rv_mean,rv_std)
% P : number of basis, r : number of multiply 
% nHr = factorial(P+r-1)/factorial(P-1)/factorial(r); 
%% Properties
[nSample,P] = size(basis);

%% Galerkin projection
E1 = zeros(P,1);
E1L = zeros(P,1);
for i = 1:P
    Pe = sum(basis(:,i))/nSample;
    PeL = sum(basis(:,i) .* (rv_mean + rv_std*XiSample))/nSample;
    E1(i,1) = Pe;
    E1L(i,1) = PeL;
end

E2 = zeros(P,P);
E2LL = zeros(P,P);
for i = 1:P
    for j = i:P
        Pe = sum(basis(:,i) .* basis(:,j))/nSample;
        PeLL = sum(basis(:,i) .* basis(:,j) .* (rv_mean + rv_std*XiSample).^2)/nSample;
        E2(i,j)= Pe; E2(j,i)= Pe;    
        E2LL(i,j)= PeLL; E2LL(j,i)= PeLL;    
    end
end
Ephi2 = diag(E2);

E3 = zeros(P,P,P);
for i = 1:P
    for j = i:P
        for l = j:P

            Pe = sum(basis(:,i) .* basis(:,j) .* basis(:,l))/nSample;

            E3(j,l,i) = Pe;
            E3(j,i,l) = Pe;
            E3(l,j,i) = Pe;
            E3(l,i,j) = Pe;
            E3(i,l,j) = Pe;
            E3(i,j,l) = Pe;
        end
    end
end
E3 = reshape(E3,P,P*P);

E4 = zeros(P,P,P,P);
for i = 1:P
    for j = i:P
        for l = j:P
            for m = l:P
                Pe = sum(basis(:,i) .* basis(:,j) .* basis(:,l) .* basis(:,m))/nSample;
                ind = perms([i,j,l,m]);
                for ii = 1:size(ind,1)
                    E4(ind(ii,1),ind(ii,2),ind(ii,3),ind(ii,4)) = Pe;
                end
            end
        end
    end
end
E4 = reshape(E4,P,P*P*P);

Galerkin.E1 = E1;
Galerkin.E1L = E1L;
Galerkin.E2 = E2;
Galerkin.E2LL = E2LL;
Galerkin.Ephi2 = Ephi2;
Galerkin.E3 = E3;
Galerkin.E4 = E4;



