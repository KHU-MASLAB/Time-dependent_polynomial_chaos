function [opts,init] = SinglePendulum(opts)
arguments % 
    %% simulation property
    opts.start = 0;
    opts.h = 1e-2; % time Step
    opts.endTime = 30;
    opts.g = 9.80665;  % m/s^2
    %% Rigid body property
    opts.L = 1; % m
    opts.radius = 0.01; % m
    opts.rho = 7.830e-6; % kg/mm^3   
    opts.th1_0 = 0;
end
%% 
opts.tspan = opts.start:opts.h:opts.endTime;
opts.A = (opts.radius^2)*pi; % m^2
opts.m1 = opts.rho*opts.A*opts.L;

%% rigid body
Rx1 = opts.L*cos(opts.th1_0)/2;    Ry1 = opts.L*sin(opts.th1_0)/2;    th1 = opts.th1_0;
opts.j1 = (1/12)*opts.m1*opts.L^2;
opts.mass = diag([opts.m1, opts.m1, opts.j1]);
%% Initialize
q = [Rx1 Ry1 th1]';
dq = zeros(3,1);
init = [q;dq];

