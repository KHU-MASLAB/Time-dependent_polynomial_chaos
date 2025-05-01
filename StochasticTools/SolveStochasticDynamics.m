function PC = SolveStochasticDynamics(Prop,init,Galerkin,BasisSample)

%% initial config. (theta1 = 0)
P = length(Galerkin.E2);
DOF = length(init)/2;
init_PC = zeros(DOF*P*2,1);
init_PC(1:P) = 0.5*Galerkin.E1L*cos(Prop.th1_0);
init_PC(P+1:2*P) = 0.5*Galerkin.E1L*sin(Prop.th1_0);
init_PC(2*P+1:3*P) = Galerkin.E1*Prop.th1_0;

%% stochastic solver
ode_opt = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t_PC,X] = ode45(@(t,x) SolveDAE_stochastic(t,x,Prop,Galerkin),Prop.tspan,init_PC,ode_opt);

%% pre process
mean_PC = X(:,1:P:DOF*2*P);
var_PC = zeros(length(t_PC),DOF*2);
for i = 1:DOF*2
    var_PC(:,i) = X(:,(i-1)*P+1:i*P).^2*Galerkin.Ephi2 - mean_PC(:,i).^2;
end

skew_PC = zeros(length(t_PC),DOF*2);
for i = 1:P
    for j = 1:P
        for k = 1:P
            a = (i-1)*P + k;
            for l = 1:DOF*2
                skew_PC(:,l) = skew_PC(:,l) + Galerkin.E3(j,a)*X(:,(l-1)*P+i).*X(:,(l-1)*P+j).*X(:,(l-1)*P+k);
            end
        end
    end
end
kurt_PC = zeros(length(t_PC),DOF*2);
for i = 2:P
    for j = 2:P
        for k = 2:P
            for m = 2:P
                a = (i-1)*P^2 + (j-1)*P + m;
                for l = 1:DOF*2
                    kurt_PC(:,l) = kurt_PC(:,l) +...
                        Galerkin.E4(k,a)*X(:,(l-1)*P+i).*X(:,(l-1)*P+j).*X(:,(l-1)*P+k).*X(:,(l-1)*P+m);
                end
            end
        end
    end
end
skew_PC = (skew_PC - 3*mean_PC.* var_PC - mean_PC.^3)./sqrt(var_PC).^3;
kurt_PC = kurt_PC./var_PC.^2;

PC.th1Sample_PC = (BasisSample * X(:,2*P+1:3*P)')';
PC.th1DotSample_PC = (BasisSample * X(:,(2+DOF)*P+1:(3+DOF)*P)')';

PC.mean_PC = mean_PC;
% PC.var_PC = sqrt(var_PC);
PC.var_PC = var_PC;
PC.skew_PC = skew_PC;
PC.kurt_PC = kurt_PC;
PC.t_PC = t_PC;


