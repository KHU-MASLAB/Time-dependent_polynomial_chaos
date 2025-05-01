function co = correlation(rve)
    [Nsample,Nvar] = size(rve);
    rv = rve - sum(rve)/Nsample;
    std_rv2 = sqrt(sum(rv.^2)/Nsample);
    co = zeros(Nvar);
    for j = 1:Nvar
        for k = 1:Nvar
            co(j,k) = (sum( rv(:,j) .* rv(:,k) )/Nsample)/std_rv2(j)/std_rv2(k);
        end
    end
end