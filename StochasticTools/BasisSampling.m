function [XiSample,Basis] = BasisSampling(nSample,D,o,dist)
% nSample : number of PC sample 
% D : number of random variable
% o : polynomial order
% dist : PDF - string
%% xi sampling
rng default
if dist == "Uniform"
    XiSample = (2 * lhsdesign(nSample,1) - 1);
elseif dist == "Normal"
    XiSample = icdf( "Normal",lhsdesign(nSample,1),0,1);
end
%% basis sampling
SingleBasis = zeros(nSample,o+1,D);
for j = 1:D
    if dist == "Uniform"
        for i = 1 : nSample
            SingleBasis(i,:,j) = Legendre_P(o,XiSample(i));
        end
    elseif dist == "Normal"
        for i = 1 : nSample
            SingleBasis(i,:,j) = hermite_F(o,XiSample(i));
        end
    end
end
%% combine
midx = Multi_Index(D,o);
P = factorial(D+o)/factorial(D)/factorial(o); % number of stochastic basis
Basis = ones(nSample, P);
for i = 1:P
    for j = 1:D
        Basis(:,i) = Basis(:,i) .* SingleBasis(:,midx(i,j)+1);
    end
end
%% Normalize
for i = 1:P
    Basis(:,i) = Basis(:,i)/sqrt(sum(Basis(:,i).^2)/nSample);
end
%% Orthogonal check
orth = zeros(P,P);
for i = 1:P
    for j = 1:P
        orth(i,j) = Basis(:,i)' * Basis(:,j) / nSample;
    end
end

end
