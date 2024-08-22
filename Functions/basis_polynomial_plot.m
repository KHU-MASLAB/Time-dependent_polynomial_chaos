function basis_polynomial_plot(xi_sample,basis_sample)

[ss,dd] = sort(xi_sample);
Nbasis = size(basis_sample,2);
close all
figure()
hold on
grid on
for i = 1:Nbasis
    plot(ss,basis_sample(dd,i),'LineWidth',2)   
end
legend('n=0','n=1','n=2','n=3','n=4')
axis([-3,3,-5,5])