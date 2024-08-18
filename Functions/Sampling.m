function [basis_sample,xi_sample] = Sampling(o,TD_sample,dist)

rng default
if dist == "Uniform"
    xi_sample = (2 * lhsdesign(TD_sample,1) - 1);
elseif dist == "Normal"
    xi_sample = icdf( "Normal",lhsdesign(TD_sample,1),0,1);
end

basis_sample = zeros(TD_sample,o+1);
if dist == "Uniform"
    for i = 1 : TD_sample
        basis_sample(i,:) = legendreP(0:o,xi_sample(i));
    end
elseif dist == "Normal"
    for i = 1 : TD_sample
        basis_sample(i,:) = hermite_F(o,xi_sample(i));
    end
end
for i = 1:o+1
    basis_sample(:,i) = basis_sample(:,i)/sqrt(sum(basis_sample(:,i).^2)/TD_sample);
end
orth = zeros(o+1,o+1);
for i = 1:o+1
    for j = 1:o+1
        orth(i,j) = basis_sample(:,i)' * basis_sample(:,j) / TD_sample;
    end
end
disp("Check orthogonality of initial basis :")
disp(orth)
end