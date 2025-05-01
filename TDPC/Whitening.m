function Basis = Whitening(phi)
[nSample,nBasis] = size(phi);
G = zeros(nBasis,nBasis);
for i = 1:nBasis
    for j = i:nBasis
        Mij = sum(phi(:,i) .* phi(:,j))/nSample;
        G(i,j) = Mij;
        G(j,i) = Mij;
    end
end
W = chol(G);
invW = W\eye(length(W));
Basis = zeros(nSample,nBasis);
for i = 1:nBasis
    for j = 1:i
        Basis(:,i) = Basis(:,i) + invW(j,i) * phi(:,j);
    end
end
