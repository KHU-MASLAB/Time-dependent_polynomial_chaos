function [phiq, gamma] = Constraint(q1,q1Dot,L)

% Rx1 = q1(1);
% Ry1 = q1(2);
th1 = q1(3);
th1Dot = q1Dot(3);

matA = [cos(th1), -sin(th1);
        sin(th1),  cos(th1)];
matB = [-sin(th1), -cos(th1);
         cos(th1), -sin(th1)];
sij = [-L/2; 0];

%% Revolute joint
phiqR = eye(2);
phiqt = matB*sij;
phiq = [phiqR,phiqt];

gamma = matA*sij*th1Dot^2;




