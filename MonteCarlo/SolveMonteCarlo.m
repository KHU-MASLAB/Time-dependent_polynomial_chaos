function out = SolveMonteCarlo(opts)

arguments
    opts.dist = "Normal";
    opts.rv_mean = 1;
    opts.rv_level = 0.05;
    opts.nSample = 10^3;
    opts.rv_name = "L";
    opts.SampleDir = '';
    opts.h = 1e-2; % time Step
    opts.endTime = 30;
    opts.th1_0 = 0;
end
opts.rv_std = opts.rv_mean*opts.rv_level;
    
%% random variable
if opts.dist == "Uniform"
    xi_M = lhsdesign(opts.nSample,1);
    rv_M = opts.rv_mean - opts.rv_std + xi_M*2*opts.rv_level;
elseif opts.dist == "Normal"
    xi_M = icdf('Normal',lhsdesign(opts.nSample,1),0,1);
    rv_M = opts.rv_mean + opts.rv_std*xi_M;
end

%% Solve dynamics system
[Prop,init] = SinglePendulum("h",opts.h,"endTime",opts.endTime,opts.rv_name,opts.rv_mean,"th1_0",opts.th1_0);
ode_opt = odeset('RelTol',1e-6,'AbsTol',1e-8);
step = length(Prop.tspan);
DOF = length(init)/2;

MC_mean = zeros(step,2*DOF);
MC_var = zeros(step,2*DOF);
MC_skew = zeros(step,2*DOF);
MC_kurt = zeros(step,2*DOF);

%% parallel works
rvName = opts.rv_name;
h = opts.h; endTime = opts.endTime; th1_0 = opts.th1_0;

tempfile = tempname(pwd); 
A = fopen(tempfile,'wt'); 

for i = 1:opts.nSample
    [Prop,init] = SinglePendulum("h",h,"endTime",endTime,rvName,rv_M(i),"th1_0",th1_0);
    [~, X] = ode45(@(t,x) SolveDAE_SinglePendulum(t,x,Prop), Prop.tspan, init, ode_opt);
    
    xi_i = [xi_M(i) zeros(1,size(X,1)-1)];
    fprintf(A,[repmat('%.6f,\t',1,size(X,1)-1) '%.6f\t' '\n'],xi_i);
    fprintf(A,[repmat('%.6f,\t',1,size(X,1)-1) '%.6f\t' '\n'],X); % Write results
    
    MC_mean = MC_mean + X;
    MC_var = MC_var + X.^2;
    MC_skew = MC_skew + X.^3;
    MC_kurt = MC_kurt + X.^4;
end
fclose(A);
tbl = readmatrix(tempfile,"Delimiter",",");
delete(tempfile)

%% save data
if isempty(opts.SampleDir)
    dirname = strrep(string(datetime),':','-');
    mkdir(dirname);
else
    dirname = opts.SampleDir;
    FileName = dirname+"\MC.txt";
    PreSim = readmatrix(FileName);
    tbl = vertcat(PreSim,tbl);
end
FileName = dirname+"\MC.txt";
% writematrix(tbl, FileName);
A = fopen(FileName,'wt'); 
fprintf(A,[repmat('%.6f,\t',1,size(tbl,2)-1) '%.6f\t' '\n'],tbl');
fclose(A);

[Prop,init] = SinglePendulum("h",opts.h,"endTime",opts.endTime,opts.rv_name,opts.rv_mean,"th1_0",opts.th1_0);
PropFileName = dirname + "\MC_Parameters.mat";


%% Calculate moments
out.t_MC = Prop.tspan;
out.mean_MC = MC_mean/opts.nSample;
out.var_MC  = MC_var/opts.nSample - out.mean_MC.^2;
out.skew_MC = (MC_skew/opts.nSample - 3*out.mean_MC.*out.var_MC - out.mean_MC.^3)./sqrt(out.var_MC).^3;
out.kurt_MC = (MC_kurt/opts.nSample - 4*out.mean_MC.*out.skew_MC/opts.nSample + 6*out.mean_MC.^2.* out.var_MC/opts.nSample...
            -3*out.mean_MC.^4)./out.var_MC.^2;

save(PropFileName,"Prop","init","opts","out");

