function Basis_combine = GenBasis(loc,D,o,dist)
% loc : quadrature points
% D : number of random variable
% o : polynomial order
% dist : PDF - string
%% basis sampling
Npoint = length(loc);
basis = zeros(Npoint,o+1,D);
for j = 1:D
    if dist == "Uniform"
        for i = 1 : Npoint
            basis(i,:,j) = Legendre_P(o,loc(i));
        end
    elseif dist == "Normal"
        for i = 1 : Npoint
            basis(i,:,j) = hermite_F(o,loc(i));
        end
    end
end
%% Combine
midx = Multi_Index(D,o);
P = factorial(D+o)/factorial(D)/factorial(o); % number of stochastic basis
Basis_combine = ones(Npoint, P,D);
for i = 1:P
    for j = 1:D
        Basis_combine(:,i,j) = basis(:,midx(i,j)+1,j);
    end
end