function M = Mass_PC(Prop,Galerkin)

P = length(Galerkin.E2);
m = Prop.mass(1,1);
j = Prop.mass(3,3);

M = zeros(3*P);
M(1:P,1:P) = m*Galerkin.E2;
M(P+1:2*P,P+1:2*P) = m*Galerkin.E2;
M(2*P+1:3*P,2*P+1:3*P) = (j/Prop.L^2)*Galerkin.E2LL;