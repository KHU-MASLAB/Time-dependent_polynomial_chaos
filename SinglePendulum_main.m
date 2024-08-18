clc; clear;
addpath('Functions');
addpath('Plot')

%% %%%%%%%%% User defined %%%%%%%%%%%%%%%

% Distribution
Distribution = "Normal";        % available option : Normal, Uniform

% model property
m = 1;                          % mass [kg]
Lb = 1;                         % mean of the pendulum length [m]
Ld = 0.05;                      % std or range of the pendulum length [m]
theta = -deg2rad(45);           % initial config.[rad]

% Monte Carlo
n_sample = 1E6;                 % Number of MC sample

% PC set
o = 4;                          % degree of polynomial basis
TD_sample = 1E5;                % Number of TD PC sample

% solver set
start = 0;                      % start time
endtime = 30;                   % end time
h = 1e-2;                       % time step

% plot flag
basis_plot = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% gauss quadrature
if Distribution == "Uniform"
    gauss_opt = '5point';
    [wts, loc] = gaussQuadrature(gauss_opt); % gauss quadrature
elseif Distribution == "Normal"
    gauss_point = 70;
    [loc,wts,~] = gengausshermquadrule2(gauss_point); % Normal dist.
    loc = loc'; wts = wts';
end
Npoint = length(wts);
%% property
g = 9.8;                % gravity (m/s)
J = (1/12)*m*Lb^2;      % moment of inertia

%% solver set
step = endtime/h;
tspan = start:h:endtime;
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
stopcrt_t = 0.1;        % Minimum time step to check update criterial

%% sampling
[basis_sample,xi_sample] = Sampling(o,TD_sample,Distribution); % [initial basis, xi sample]
D = 1;          % dimension of random variable
P = factorial(D+o)/factorial(D)/factorial(o); % number of stochastic basis
if basis_plot ==1
    basis_polynomial_plot(xi_sample,basis_sample)
end
filename = "SP_L"+num2str(Lb)+"_dL" +num2str(Ld)+"_P"+num2str(o)+"_MC"...
            + "_PC"+num2str(TD_sample,'%.0E')+"_t"+num2str(endtime)+"_h"+num2str(h);
filenameMC = "SP_L"+num2str(Lb)+"_dL" +num2str(Ld)+"_MC"+num2str(n_sample,'%.0E')...
    + "_t"+num2str(endtime)+"_h"+num2str(h);
%% %%%%%%%%%%%%%%%%%%% Monte Carlo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist("result/MC_"+filenameMC+".mat",'file'))
    load("result/MC_"+filenameMC+".mat");
    fprintf("MC Imported /MC_"+filenameMC+"\n")
else
    tic
    [t_MC,mean_MC,var_MC,skew_MC,kurt_MC,thSample_MC,dthSample_MC]...
    =  MC(Distribution,g,m,Lb,Ld,h,n_sample,endtime,tspan,theta,opts);
    runtime_MC = toc;
    fprintf('MC finished, runtime : %f\n',runtime_MC)
    save("result/MC_"+filenameMC,"t_MC","mean_MC","var_MC","skew_MC","kurt_MC",...
        "thSample_MC","dthSample_MC");
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%  PC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist("result/PC_"+filename+".mat",'file'))
    load("result/PC_"+filename+".mat");
    fprintf("PC imported /PC_"+filename+"\n")
else
    tic
    [E, Ephi2, E_L2, E1, E2, E3, E4,t_PC,mean_PC,var_PC,skew_PC,kurt_PC,XSample_PC,YSample_PC,thSample_PC,dthSample_PC]...
    = PC(Distribution,Npoint,P,o,h,loc,wts,Lb,Ld,theta,tspan,opts,m,g,basis_sample);
    runtime_PC = toc;
    fprintf('PC finished, runtime : %f\n',runtime_PC)
    save("result/PC_"+filename,"t_PC","mean_PC","var_PC","skew_PC","kurt_PC",...
        "thSample_PC","XSample_PC","YSample_PC","dthSample_PC");
end
%% %%%%%%%%%%%%%%%%%%%%%%%%% TD PC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist("result/TDPC_"+filename+".mat",'file'))
    load("result/TDPC_"+filename+".mat");
    fprintf("TD-PC Imported/TDPC_"+filename+"\n")
else
    tic
    [t_TDPC,mean_TDPC,var_TDPC,skew_TDPC,kurt_TDPC,thSample_TDPC,XSample_TDPC,YSample_TDPC,dthSample_TDPC,Nupdate]...
     = TDPC(endtime,stopcrt_t,tspan,Ephi2,E,E1,E2,E3,E4,E_L2,basis_sample,xi_sample,...
            theta,h,o,Lb,Ld,g,m,TD_sample);
    runtime_TDPC = toc;
    fprintf('TD-PC finished, runtime : %f, Nupdate : %f\n',runtime_TDPC,Nupdate)
    save("result/TDPC_"+filename,"t_TDPC","mean_TDPC","var_TDPC","skew_TDPC","kurt_TDPC",...
        "thSample_TDPC","XSample_TDPC","YSample_TDPC","dthSample_TDPC");
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% result plot
% plotall(t_MC,t_PC,t_TDPC,mean_MC,mean_PC,mean_TDPC,"mean")
% plotall(t_MC,t_PC,t_TDPC,var_MC,var_PC,var_TDPC,"Variance")
plotall(t_MC,t_PC,t_TDPC,real(skew_MC),skew_PC,skew_TDPC,"Skewness")
% plotall(t_MC,t_PC,t_TDPC,real(kurt_MC),kurt_PC,kurt_TDPC,"Kurtosis")
%% Error plot
% plotall_Error(t_PC,t_TDPC,mean_MC,mean_PC,mean_TDPC,"mean",Distribution,filename)
% plotall_Error(t_PC,t_TDPC,var_MC,var_PC,var_TDPC,"Variance",Distribution,filename)
plotall_Error(t_PC,t_TDPC,real(skew_MC),skew_PC,skew_TDPC,"Skewness",Distribution,filename)
% plotall_Error(t_PC,t_TDPC,real(kurt_MC),kurt_PC,kurt_TDPC,"Kurtosis",Distribution,filename)
%% Distribution
close all
tt = 20;
plot_distribution(t_PC,t_TDPC,thSample_MC,thSample_PC,thSample_TDPC,tt,"$\theta$")
plot_distribution(t_PC,t_TDPC,dthSample_MC,dthSample_PC,dthSample_TDPC,tt,"d$\theta$")
