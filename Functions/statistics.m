function [mean_temp,var_temp,skew_temp,kurt_temp,rvth,rvdth] =...
    statistics(TDbasis,DOF,ye,P,lent)
Nsample = size(TDbasis,1);
mu_basis = sum(TDbasis)/Nsample;
mean_temp = zeros(2*DOF,lent);
var_temp  = zeros(2*DOF,lent);
skew_temp = zeros(2*DOF,lent);
kurt_temp = zeros(2*DOF,lent);
for j = 1:2*DOF
    mean_temp(j,:) =  mu_basis * ye(1:lent,1+P*(j-1):j*P)';
    rv_temp = TDbasis * ye(1:lent,1+P*(j-1):j*P)';
    mu2 = sum(rv_temp.^2)/Nsample;
    mu3 = sum(rv_temp.^3)/Nsample;
    mu4 = sum(rv_temp.^4)/Nsample;
    var_temp(j,:) = mu2 - mean_temp(j,:).^2;
    skew_temp(j,:) = real((mu3 - 3*mean_temp(j,:).*var_temp(j,:) - mean_temp(j,:).^3)./(sqrt(var_temp(j,:)).^3));
    kurt_temp(j,:) = real((mu4 - 4*mean_temp(j,:).*mu3 + 6*mean_temp(j,:).^2.*mu2 - 3*mean_temp(j,:).^4)./var_temp(j,:).^2);
end

rvth  =  TDbasis * ye(1:end-2,1+2*P:3*P)';
rvdth  =  TDbasis * ye(1:end-2,1+5*P:6*P)';
        
end