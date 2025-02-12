function [Galerkin,NewBasis,init_TDPC] = updategPC_GramSchmidt_sampling(opts,rv,Xe,XiSample,BasisSample)

%% property
[nSample,P] = size(BasisSample);
Nvar = size(rv,2);
order = opts.o;
nBasis = factorial(Nvar + order)/factorial(Nvar)/factorial(order);
DOF = length(Xe)/P;

%% sample variables
X = reshape(Xe,P,DOF);
rvX = BasisSample*X;
rv = rv - sum(rv)/nSample;

%% multi index
midx = Multi_Index(Nvar,order);

%% multivariate monomials
Xbeta = MultiMonomial(rv,midx);

%% modified Gram-Schmidt 
% phi = ModifiedGS(Xbeta);
phi = GS(Xbeta);

%% Whitening transformation
NewBasis = Whitening(phi);

%% Galerkin projection
Galerkin = GalerkinUpdate(NewBasis,XiSample,opts.rv_mean,opts.rv_std);

%% update Initial Condition

init_TDPC = zeros(nBasis,DOF);
for j = 1:nBasis
    init_TDPC(j,:) = sum(rvX .* NewBasis(:,j))/nSample;    
end
init_TDPC = reshape(init_TDPC,nBasis*DOF,1);

end