function phi = ModifiedGS(Xbeta)

[nSample,nBasis] = size(Xbeta);
phi = zeros(nSample,nBasis);

phi(:,1) = Xbeta(:,1) / sqrt(sum(Xbeta(:,1))/nSample);
for i = 2:nBasis
    phi(:,i) = Xbeta(:,i);
    for j = 1:i-1
        phi(:,i) = phi(:,i) - ((phi(:,j)'*phi(:,i))/nSample) * phi(:,j);
    end
    phi(:,i) = phi(:,i) / sqrt(sum(phi(:,i).^2)/nSample);
end