function [t_MC,mean_MC,var_MC,skew_MC,Kurt_MC,thSample_MC,dthSample_MC] =...
    MC(dist,g,m,Lb,Ld,h,n_sample,endtime,tspan,theta,opts)
sol = struct();
if dist == "Uniform"
    L_M = Lb - Ld + lhsdesign(n_sample,1)*2*Ld;
elseif dist == "Normal"
    L_M = Lb + Ld*icdf('Normal',lhsdesign(n_sample,1),0,1);
end
sol.MC_mean = zeros(length(tspan),6);
sol.MC_var = zeros(length(tspan),6);
sol.MC_skew = zeros(length(tspan),6);
sol.MC_kurt = zeros(length(tspan),6);
sol.th = zeros(endtime+1,n_sample);
sol.dth = zeros(endtime+1,n_sample);
for i = 1:n_sample
    init_MC = [L_M(i)*cos(theta); L_M(i)*sin(theta); theta; 0;0;0];
    [t_MC,X_MC] = ode45(@(t,x) SP_ode_MC(t,x,L_M(i),g,m),tspan,init_MC,opts);

    sol.MC_mean = sol.MC_mean + X_MC;
    sol.MC_var = sol.MC_var + X_MC.^2;
    sol.MC_skew = sol.MC_skew + X_MC.^3;
    sol.MC_kurt = sol.MC_kurt + X_MC.^4;
    sol.th(:,i) = X_MC(1:1/h:length(tspan),3);
    sol.dth(:,i) = X_MC(1:1/h:length(tspan),6);
end

mean_MC = sol.MC_mean/n_sample;
var_MC  = sol.MC_var/n_sample - mean_MC.^2;
skew_MC = (sol.MC_skew/n_sample - 3*mean_MC.*var_MC - mean_MC.^3)./sqrt(var_MC).^3;
% skew_MC = (sol.MC_skew/n_sample - 3*mean_MC.*sol.MC_var/n_sample + 2*mean_MC.^3);
Kurt_MC = (sol.MC_kurt/n_sample - 4*mean_MC.*sol.MC_skew/n_sample + 6*mean_MC.^2.*sol.MC_var/n_sample...
            -3*mean_MC.^4)./var_MC.^2;
% Kurt_MC = (sol.MC_kurt/n_sample - 4*mean_MC.*sol.MC_skew/n_sample + 6*mean_MC.^2*sol.MC_var/n_sample...
%             -3*mean_MC.^4);
thSample_MC = sol.th';
dthSample_MC = sol.dth';





