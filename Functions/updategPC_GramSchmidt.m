function [new_basis,Ephi2,E,E1, E2, E3, E_L2,init_TDPC] = updategPC_GramSchmidt(t,rv,Xe,xi_sample,o,Lb,Ld,TDbasis)


Nsample = length(TDbasis);
P = size(TDbasis,2);
x   =  Xe(1:P);
y   =  Xe(P+1:2*P);
th  =  Xe(2*P+1:3*P);
dx   =  Xe(3*P+1:4*P);
dy   =  Xe(4*P+1:5*P);
dth  =  Xe(5*P+1:6*P);

% rv = [ TDbasis*x', TDbasis*dx',TDbasis*y',TDbasis*dy'];
% rv = [ TDbasis*th'];

mean = sum(rv)/Nsample;
% varMC = sum(rv.^2)/Nsample - (mean).^2;
% varPC = [sum(th(2:end).^2),sum(dth(2:end).^2)];
rv = rv - mean;

rvX = [TDbasis*x', TDbasis*y', TDbasis*th',TDbasis*dx', TDbasis*dy', TDbasis*dth'];
mean_X = sum(rvX)/Nsample;
varMC = sum(rvX.^2)/Nsample - (mean_X).^2;
varPC = [sum(x(2:end).^2),sum(y(2:end).^2),sum(th(2:end).^2),sum(dx(2:end).^2),sum(dy(2:end).^2),sum(dth(2:end).^2)];
% rvX = rvX - mean_X;

% 
% close all
% figure();
% hold on
% for i = 2:size(rv,2)
% subplot(3,2,i-1)
% plot(rv(:,1),rv(:,i),'o')
% end
% 
% figure();
% hold on
% for i = 2:size(rv,2)
% subplot(3,2,i-1)
% histogram(rv(:,i))
% end

%% property
D = 0;
Nvar = size(rv,2);
order = o;
Nbasis = factorial(Nvar + order)/factorial(Nvar)/factorial(order);


%% multi index
midx = zeros(Nbasis,Nvar);
a = 1;
for i = 0:order
    if i == 0
        midx(1,:) = zeros(1,Nvar);
    else
        num = (factorial(Nvar + i-1)/factorial(Nvar)/factorial(i-1))*(Nvar/i);
        midx(a+1:a+num,:) = multi_idx(i,Nvar);
        a = a + num;
    end
end
% gro = [];
% for i = 5:Nbasis
%     if sum(midx(i,:)) == midx(i,1)
%         gro = [gro,i];
%     end
% end
% midx(gro',:) = [];
% midx(end,:) = [];
% Nbasis = size(midx,1);

%% Gram-Schmidt
Xbeta = zeros(Nsample, Nbasis);
for i = 1:Nbasis
%     Pe = xi_sample.^midx(i,1);
    Pe = ones(Nsample,1);
    for k = 1:Nvar
        Pee = rv(:,k).^midx(i,k);
        Pe = Pe .* Pee;
    end
    Xbeta(:,i) = Pe;
end

std_X = sqrt(sum(Xbeta.^2)/Nsample);
Corr = zeros(Nbasis,Nbasis);
for i = 1: Nbasis
    for j = 1:Nbasis
        Corr(i,j) = (sum( Xbeta(:,i) .* Xbeta(:,j) )/Nsample)/std_X(i)/std_X(j);
    end
end

%%
new_basis = zeros(Nsample,Nbasis);
new_basis(:,1) = 1;
U_sample = zeros(Nbasis,Nbasis);
L_sample = zeros(Nbasis,Nbasis);
for i = 2:Nbasis
    new_basis(:,i) = Xbeta(:,i);
    for j = 1:i-1
        U = sum(Xbeta(:,i) .* new_basis(:,j))/Nsample;
        L = sum(new_basis(:,j).^2 )/Nsample;
        U_sample(i,j) = U;
        L_sample(i,j) = L;
        new_basis(:,i) = new_basis(:,i) - new_basis(:,j)*(U/L);
    end
end
% Normalize
for i = 1:Nbasis
    new_basis(:,i) = new_basis(:,i) / sqrt(sum(new_basis(:,i) .* new_basis(:,i))/Nsample);
end
%%
new_basis = zeros(Nsample,Nbasis);

new_basis(:,1) = Xbeta(:,1) / sqrt(sum(Xbeta(:,1))/Nsample);
for i = 2:Nbasis
    new_basis(:,i) = Xbeta(:,i);
    for j = 1:i-1
        new_basis(:,i) = new_basis(:,i) - ((new_basis(:,j)'*new_basis(:,i))/Nsample) * new_basis(:,j);
    end
    new_basis(:,i) = new_basis(:,i) / sqrt(sum(new_basis(:,i).^2)/Nsample);
end

%% Gram matrix
M = zeros(Nbasis,Nbasis);
for i = 1:Nbasis
    for j = i:Nbasis
        Mij = sum(Xbeta(:,i) .* Xbeta(:,j))/Nsample;
        M(i,j) = Mij;
        M(j,i) = Mij;
    end
end
d = eig(M);
rk = rank(M);
% if rk<Nbasis
%     [U,S,V] = svd(M);
%     % S = S + (d(Nbasis-rk)+1e-5)*eye(Nbasis);
%     % Mmodi = U * S * V';
%     % [U,S,V] = svd(Mmodi);
%     Mmodi = M + (d(Nbasis-rk)+1e-5)*eye(Nbasis);
%     [U,S,V] = svd(Mmodi);
% %     rank(Mmodi);
% %     eig(Mmodi);
%     R = chol(Mmodi);
%     invR = R\eye(length(R));
%     psi = zeros(Nsample,Nbasis);
%     for i = 1:Nbasis
%         for j = 1:i
%             psi(:,i) = psi(:,i) + invR(j,i) * Xbeta(:,j);
%         end
%     end
%     Epsi = zeros(Nbasis,Nbasis);
%     for i = 1:Nbasis
%         for j = i:Nbasis
%             Epsi(i,j) = sum(psi(:,i) .* psi(:,j))/Nsample;
%             Epsi(j,i) = Epsi(i,j);
%         end
%     end
% end
% figure;
% plot(1:Nbasis,d)
% rk = 4;
% R = chol(M);
% invR = R\eye(length(R));
% psi = zeros(Nsample,Nbasis);
% for i = 1:Nbasis
%     for j = 1:i
%         psi(:,i) = psi(:,i) + invR(j,i) * Xbeta(:,j);
%     end
% end
% new_basis = new_basis(:,1:rk);
% Nbasis = rk;
%% eliminate non orthogonal basis
orth = zeros(Nbasis,Nbasis);
for i = 1:Nbasis
    for j = i:Nbasis
        orth(i,j) = sum(new_basis(:,i) .* new_basis(:,j))/Nsample;
        orth(j,i) = orth(i,j);
    end
end
ark = Nbasis;
for i = 2:Nbasis
        Ec = orth(1:i,1:i);
        if norm(Ec) > 1 + 1e-4
            ark = i-1;
            new_basis = new_basis(:,1:ark);
            Nbasis = size(new_basis,2);
            break;
        end        
end
%% plot basis
% close all
% figure();
% hold on
% grid on
% for i = 1:ark
%     plot(rv(:,1), new_basis(:,i),'*')
% end
% axis([-1,1,-3,3])
%% Galerkin projection
E = zeros(Nbasis,Nbasis);
for i = 1:Nbasis
    for j = i:Nbasis
        E(i,j) = sum(new_basis(:,i) .* new_basis(:,j))/Nsample;
        E(j,i) = E(i,j);
    end
end
Ephi2 = diag(E);

E_L2 = zeros(Nbasis,Nbasis);
for i = 1:Nbasis
    for j = 1:Nbasis
        Mij = sum(new_basis(:,i) .* new_basis(:,j) .* (Lb + Ld * xi_sample).^2)/Nsample;
        E_L2(i,j) = Mij;
        E_L2(j,i) = Mij;
    end
end

E3 = zeros(Nbasis,1);
for i = 1:Nbasis
    Pe = sum(new_basis(:,i))/Nsample;
    E3(i) = Pe;
end

E1 = zeros(Nbasis,Nbasis^2);
for i = 1:Nbasis
    E1e = zeros(Nbasis,Nbasis);
    for j = 1:Nbasis
        for l = j:Nbasis
            Mij = sum(new_basis(:,i) .* new_basis(:,j) .* new_basis(:,l))/Nsample;
            E1e(j,l) = Mij;
            E1e(l,j) = Mij;
         end
    end
    E1(:,Nbasis*(i-1)+1 : Nbasis*i) = E1e;
end

E2 = zeros(Nbasis,Nbasis^3);
for i = 1:Nbasis
    for ii = 1:Nbasis
        ME2e = zeros(Nbasis,Nbasis);
        for j = 1:Nbasis
            for l = j:Nbasis
                Mij = sum(new_basis(:,i) .* new_basis(:,ii) .* new_basis(:,j) .* new_basis(:,l))/Nsample;
                ME2e(j,l) = Mij;
                ME2e(l,j) = Mij;
             end
        end
        n1 = Nbasis^2*(i-1) + Nbasis*(ii-1) + 1;
        n2 = Nbasis^2*(i-1) + Nbasis*(ii);
        E2(:,n1:n2) = ME2e;
    end
end

%% update Initial Condition
% M_init = zeros(Nbasis,6);
% init = zeros(Nbasis,6);
% for i = 1:Nbasis
%     Pe = (loc.^midx(i,1))' * wts * W;
%     Pe2 = (loc.^midx(i,1) .* (Lb + Ld * loc))' * wts * W;
%     
%     rv1 = sum(rv(:,1).^midx(i,D+1))/Nsample;
%     rv2 = sum(rv(:,2).^midx(i,D+2))/Nsample;
% %     rv3 = sum(rv(:,3).^midx(i,D+3))/Nsample;
% %     rv4 = sum(rv(:,4).^midx(i,D+4))/Nsample;
%     %rve = rv1 * rv2 * rv3 * rv4;
% %     rve = rv1 * rv2;
%     
% %     Pe_X = rve * sum(rvX(:,1))/Nsample;
% %     Pe_Y = rve * sum(rvX(:,2))/Nsample;
% %     Pe_th = rve * sum(rvX(:,3))/Nsample;
% %     Pe_dX = rve * sum(rvX(:,4))/Nsample;
% %     Pe_dY = rve * sum(rvX(:,5))/Nsample;
% %     Pe_dth = rve * sum(rvX(:,6))/Nsample;
% 
% %      Pe_X = sum(rv(:,1).^midx(i,2) .* rvX(:,1))/Nsample * rv2;
% %      Pe_Y = sum(rv(:,1).^midx(i,2) .* rvX(:,2))/Nsample * rv2;
% %      Pe_dX = sum(rv(:,2).^midx(i,3) .* rvX(:,4))/Nsample * rv1;
% %      Pe_dY = sum(rv(:,2).^midx(i,3) .* rvX(:,5))/Nsample * rv1;
%      
%     Pe_X = Pe2 * sum(rv(:,1).^midx(i,D+1) .* cos(rvX(:,3)))/Nsample * rv2;
%     Pe_Y = Pe2 * sum(rv(:,1).^midx(i,D+1) .* sin(rvX(:,3)))/Nsample * rv2;
%     
%     Pe_dX =  Pe2 * sum(rv(:,1).^midx(i,D+1) .* (-sin(rvX(:,3))) )/Nsample * sum(rv(:,2).^midx(i,D+2) .* rvX(:,6))/Nsample;
%     Pe_dY =  Pe2 * sum(rv(:,1).^midx(i,D+1) .* cos(rvX(:,3)) )/Nsample * sum(rv(:,2).^midx(i,D+2) .* rvX(:,6))/Nsample;
% 
%     Pe_th  = Pe * sum(rv(:,1).^midx(i,D+1) .* rvX(:,3))/Nsample * rv2;
%     Pe_dth = Pe * sum(rv(:,2).^midx(i,D+2) .* rvX(:,6))/Nsample * rv1;
%     
%     M_init(i,:) =  [Pe_X, Pe_Y, Pe_th, Pe_dX, Pe_dY, Pe_dth];
% end
% for i = 1:6
%     init(:,i) = sum(invR .* M_init(:,i))';
% end
% 
% init_TDPC = zeros(Nbasis,6);
% 
% for j = 1:Nbasis
% %     init_TDPC(j,:) = mean_X * E3(j) + init(j,:);
%     init_TDPC(j,:) = mean_X * 0 + init(j,:);
% end
% init_TDPC = reshape(init_TDPC,Nbasis*6,1);
init_TDPC = zeros(Nbasis,6);
% for i = 1:6
%     for j = 1:Nbasis
%         init_TDPC(j,i) = sum(rvX(:,i) .* new_basis(:,j))/Nsample;
%     end
% end
for j = 1:Nbasis
    
    init_TDPC(j,1) = sum( (Lb + Ld * xi_sample) .* cos(rvX(:,3)) .* new_basis(:,j))/Nsample;
    init_TDPC(j,2) = sum( (Lb + Ld * xi_sample) .* sin(rvX(:,3)) .* new_basis(:,j))/Nsample;
    init_TDPC(j,3) = sum(rvX(:,3) .* new_basis(:,j))/Nsample;
    
    init_TDPC(j,4) = sum( -(Lb + Ld * xi_sample) .* sin(rvX(:,3)) .* rvX(:,6) .* new_basis(:,j))/Nsample;
    init_TDPC(j,5) = sum(  (Lb + Ld * xi_sample) .* cos(rvX(:,3)) .* rvX(:,6) .* new_basis(:,j))/Nsample;
    init_TDPC(j,6) = sum(rvX(:,6) .* new_basis(:,j))/Nsample;
    
end

init_TDPC = reshape(init_TDPC,Nbasis*6,1);

end