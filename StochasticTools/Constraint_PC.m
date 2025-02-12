function [phiq, gamma] = Constraint_PC(q1,q1Dot,Galerkin)

P = length(Galerkin.E2);
Rx1 = q1(1:P);
Ry1 = q1(P+1:2*P);
% th1 = q1(3);
th1Dot = q1Dot(2*P+1:3*P);

%% Revolute joint
phiqR = kron(eye(2),Galerkin.E2);
phiqt = [multiply(Galerkin.E3,Ry1,1); 
        -multiply(Galerkin.E3,Rx1,1)];
phiq = [phiqR,phiqt];

th1Dotth1Dot = kron(th1Dot,th1Dot);

gamma = [-Galerkin.E4 * kron(Rx1,th1Dotth1Dot);
         -Galerkin.E4 * kron(Ry1,th1Dotth1Dot)];





