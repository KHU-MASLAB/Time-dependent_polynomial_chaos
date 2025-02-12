%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Stochastic multibody dynamics using time-dependent polynomial chaos 
%      Feb. 08, 2025
%
%      Seok-Hee Han
%      Modeling and Simulation (M&S) lab
%      Department of Mechanical Engineering, Kyung Hee University
%      ygj03020@khu.ac.kr
%
%      Single pendulum problem
%      Uncertainty quantification, multibody dynamics, 
%      time dependent polynomial chaos, long time evaluation
%
%      Reference:
%
%      S.H. Han, H.S. Choi, J.G. Kim, Updating polynomial chaos basis for addressing 
%      long time evaluation of multibody systems with uncertainties
%   
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
restoredefaultpath;
addpath('Functions');
addpath('MonteCarlo');
addpath('PlotTools');
addpath('StochasticTools');
addpath('TDPC');
set(0,'DefaultFigureWindowStyle','docked')
%% Deterministic model
[Prop,init] = SinglePendulum();
Det = SolveDeterministic(Prop,init);
video(Det.t,Det.disp,Prop)

%% MonteCarlo
% MC = SolveMonteCarlo_par("rv_mean",2,"rv_level",0.05,"nSample",500,"nWorkers",8,"th1_0",-deg2rad(45));
% MC = SolveMonteCarlo("rv_mean",2,"rv_level",0.05,"nSample",500,"th1_0",-deg2rad(45),"endTime",30);
MC = ReadMonteCarlo("MC_Sample500_L2_dL0.05_th45");

%% Polynomial chaos
PC = SolvePC("rv_level",0.05,"o",4,"th1_0",-deg2rad(45),"rv_mean",2);

%% Time dependent polynomial chaos
TDPC = SolvePC_TDPC("rv_mean",2,"rv_level",0.05,"th1_0",-deg2rad(45),"o",4,"nSample",10000,"stopcrt_t",0.1);

%% result plot
% plot_moments(MC,PC,TDPC,"mean")
plot_moments(MC,PC,TDPC,"var")
% plot_moments(MC,PC,TDPC,"skew")
% plot_moments(MC,PC,TDPC,"kurt")

%% Error plot
plotall_Error(MC,PC,TDPC,"mean")
% plotall_Error(MC,PC,TDPC,"var")
% plotall_Error(MC,PC,TDPC,"skew")
% plotall_Error(MC,PC,TDPC,"kurt")
