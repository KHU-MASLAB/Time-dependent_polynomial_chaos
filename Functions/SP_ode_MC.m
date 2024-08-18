function [dxdt,Y] = SP_ode_MC(t,x,L,g,m)
dxdt=zeros(6,1);
X=x(1:3);
dX=x(4:6);
dxdt(1:3)=dX;
J = (1/12)*m*L^2;      % moment of inertia
% A=[m 0 0 1 0;
%      0 m 0 0 1;
%      0 0 J Lb*sin(X(3)) -Lb*cos(X(3));
%      1 0 Lb*sin(X(3)) 0 0;
%      0 1 -Lb*cos(X(3))  0 0];
% 
% B=[0; -m*g; 0; -Lb*cos(X(3))*dX(3)^2; -Lb*sin(X(3))*dX(3)^2];
A=[ m 0 0 1 0;
    0 m 0 0 1;
    0 0 J X(2) -X(1);
    1 0 X(2) 0 0;
    0 1 -X(1)  0 0];
% B=[0; -m*g; 0;-X(1)*(dX(1)^2+dX(2)^2)/L^2; -X(2)*(dX(1)^2+dX(2)^2)/L^2];
B = [0; -m*g; 0; -X(1)*dX(3)^2;-X(2)*dX(3)^2];
% 
% A=[m 0 0 1 0;
%     0 m 0 0 1;
%     0 0 J X(2) -X(1);
%     1 0 X(2) 0 0;
%     0 1 -X(1)  0 0];
% B=[0; -m*g; 0;-X(1)*(dX(1)^2+dX(2)^2)/Lb^2; -X(2)*(dX(1)^2+dX(2)^2)/Lb^2];

Y=linsolve(A,B);
dxdt(4:6)=Y(1:3);

