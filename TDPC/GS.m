function phi = GS(Xbeta)

[nSample,nBasis] = size(Xbeta);
phi = zeros(nSample,nBasis);

phi(:,1) = 1;
U_sample = zeros(nBasis,nBasis);
L_sample = zeros(nBasis,nBasis);
for i = 2:nBasis
    phi(:,i) = Xbeta(:,i);
    for j = 1:i-1
        U = sum(Xbeta(:,i) .* phi(:,j))/nSample;
        L = sum(phi(:,j).^2 )/nSample;
        U_sample(i,j) = U;
        L_sample(i,j) = L;
        phi(:,i) = phi(:,i) - phi(:,j)*(U/L);
    end
end
% Normalize
for i = 1:nBasis
    phi(:,i) = phi(:,i) / sqrt(sum(phi(:,i) .* phi(:,i))/nSample);
end