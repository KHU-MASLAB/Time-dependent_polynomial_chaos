function out = SolveDeterministic(Prop,init)

ode_opt = odeset('RelTol',1e-6,'AbsTol',1e-8);

DOF = length(init)/2;
nConst = 2; % number of constraint

[out.t, X] = ode45(@(t,x) SolveDAE_SinglePendulum(t,x,Prop), Prop.tspan, init, ode_opt);

out.disp = X(:,1:DOF);
out.vel = X(:,DOF+1:2*DOF);
out.acc = zeros(length(Prop.tspan),DOF);
out.lambda = zeros(length(Prop.tspan),nConst);

for i = 1:length(out.t)
    [~,Y_temp] = SolveDAE_SinglePendulum(out.t(i),X(i,:)',Prop);
    out.acc_Direct(i,:) = Y_temp(1:DOF)';
    out.lambda_Direct(i,:) = Y_temp(DOF+1:end)';  
end




