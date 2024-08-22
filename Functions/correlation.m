function [rve,co] = correlation(rve)
    [Nsample,Nvar] = size(rve);
    rve = rve - sum(rve)/Nsample;
    std_rv2 = sqrt(sum(rve.^2)/Nsample);
    co = zeros(size(rve,2));
    for j = 1:Nvar
        for k = 1:Nvar
            co(j,k) = (sum( rve(:,j) .* rve(:,k) )/Nsample)/std_rv2(j)/std_rv2(k);
        end
    end
end