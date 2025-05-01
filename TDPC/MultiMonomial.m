function Xbeta = MultiMonomial(rv,midx)

[nSample,nVar] = size(rv);
nBasis = size(midx,1);

Xbeta = ones(nSample, nBasis);
for i = 1:nBasis
    for j = 1:nVar
        Xbeta(:,i) = Xbeta(:,i) .* rv(:,j).^midx(i,j);
    end
end

G = zeros(nBasis,nBasis);
for i = 1:nBasis
    for j = 1:nBasis
        Mij = sum(Xbeta(:,i) .* Xbeta(:,j))/nSample;
        G(i,j) = Mij;
    end
end
fprintf("condition number of Gram matrix : %d\n",cond(G))
