function [E, Ephi2, E_L2, E1, E2, E3, E4,t_PC,mean_PC,var_PC,skew_PC,kurt_PC,...
    thSample_PC,dthSample_PC] = ...
    PC(dist,Npoint,P,o,h,loc,wts,Lb,Ld,theta,tspan,opts,m,g,basis_sample)

% Galerkin projection
basis = zeros(Npoint,P);
if dist == "Uniform"
    for i = 1 : Npoint
     basis(i,:) = legendreP(0:o,loc(i));
    end
elseif dist == "Normal"
    for i = 1 : Npoint
         basis(i,:) = hermite_F(o,loc(i));
    end
end
[E, Ephi2, E_L2, E1, E2, E3, E4] = PC_Pre(dist,basis,loc,wts,Lb,Ld);
DOF = 3;

% RK4
x_i = [ E4*cos(theta); E4*sin(theta); theta; zeros(P-1,1)];
xD_i = zeros(DOF*P,1);
init_PC = [x_i;xD_i];

[t_PC,X] = ode45(@(t,x) SP_ode_PC(t,x,g, m, P, E, E1, E2, E3, E_L2, DOF),tspan,init_PC,opts);

mean_PC = X(:,1:P:6*P);
var_PC = zeros(length(t_PC),6);
for i = 1:6
    var_PC(:,i) = X(:,(i-1)*P+1:i*P).^2*Ephi2 - mean_PC(:,i).^2;
end

skew_PC = zeros(length(t_PC),6);
for i = 1:P
    for j = 1:P
        for k = 1:P
            a = (i-1)*P + k;
            for l = 1:6
                skew_PC(:,l) = skew_PC(:,l) + E1(j,a)*X(:,(l-1)*P+i).*X(:,(l-1)*P+j).*X(:,(l-1)*P+k);
            end
        end
    end
end
kurt_PC = zeros(length(t_PC),6);
for i = 2:P
    for j = 2:P
        for k = 2:P
            for m = 2:P
                a = (i-1)*P^2 + (j-1)*P + m;
                for l = 1:6
                    kurt_PC(:,l) = kurt_PC(:,l) +...
                        E2(k,a)*X(:,(l-1)*P+i).*X(:,(l-1)*P+j).*X(:,(l-1)*P+k).*X(:,(l-1)*P+m);
                end
            end
        end
    end
end
skew_PC = (skew_PC - 3*mean_PC.* var_PC - mean_PC.^3)./sqrt(var_PC).^3;
kurt_PC = kurt_PC./var_PC.^2;

thSample_PC = basis_sample * X(1:1/h:end,2*P+1:3*P)';
dthSample_PC = basis_sample * X(1:1/h:end,5*P+1:6*P)';



