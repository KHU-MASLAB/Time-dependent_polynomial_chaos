function out = ReadMonteCarlo(FoldName)

%% Load data
FileName = FoldName + "\MC.txt";
ParaFileName = FoldName + "\MC_Parameters.mat";
MC = readmatrix(FileName,"Delimiter",",");
load(ParaFileName);
nSample = size(MC,1)/(length(init) + 1);
fprintf("Number of samples : %d \n",nSample)

%% Calculate moments
step = length(Prop.tspan);
DOF = length(init);
MC_mean = zeros(step,DOF);
MC_var = zeros(step,DOF);
MC_skew = zeros(step,DOF);
MC_kurt = zeros(step,DOF);
out.th1Sample_MC = zeros(step,nSample);
out.th1DotSample_MC = zeros(step,nSample);
out.xi = zeros(1,nSample);

% figure(1);
% hold on
for i = 1:nSample
    n1 = (i-1)*(DOF+1) + 1;
    n2 = i*(DOF+1);
    out.xi(i) = MC(n1,1);
    nMC = MC(n1+1:n2,:)';
    MC_mean = MC_mean + nMC;
    MC_var = MC_var + nMC.^2;
    MC_skew = MC_skew + nMC.^3;
    MC_kurt = MC_kurt + nMC.^4;
    out.th1Sample_MC(:,i) = nMC(:,3);
    out.th1DotSample_MC(:,i) = nMC(:,3+DOF/2);    
%     plot(Prop.tspan,nMC(:,11))
end

out.t_MC = Prop.tspan;
out.mean_MC = MC_mean/nSample;
% out.var_MC  = sqrt(MC_var/nSamples - out.mean_MC.^2);
out.var_MC  = MC_var/nSample - out.mean_MC.^2;
out.skew_MC = (MC_skew/nSample - 3*out.mean_MC.*out.var_MC - out.mean_MC.^3)./sqrt(out.var_MC).^3;
out.kurt_MC = (MC_kurt/nSample - 4*out.mean_MC.*out.skew_MC/nSample + 6*out.mean_MC.^2.* out.var_MC/nSample...
            -3*out.mean_MC.^4)./out.var_MC.^2;

