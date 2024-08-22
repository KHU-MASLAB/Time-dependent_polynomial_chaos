function [dxdt,Y] = SP_ode_TDPC(t,X,g, m, E, E1, E2, E3, E_L2, DOF)

P = length(E);
dxdt = zeros(2*DOF*P,1);
dxdt(1:DOF*P) = X(DOF*P+1:2*DOF*P);

% A1=[];A2=[];
A3=[];A4=[];
B1=[];B2=[];B3=[];B4=[];
D1=[]; D2=[]; D3=[];
for i=1:P
%     A1=[A1;X(i)*X(P+1:2*P)];        % x y
%     A2=[A2;X(i)*X(1:P)];            % x x

    A3 = [A3;X(P*DOF+i)*X(P*DOF+1:P*(DOF+1))];          % Vx Vx
    A4 = [A4;X(P*(DOF+1)+i)*X(P*(DOF+1)+1:P*(DOF+2))];  % Vy Vy
    D1 = [D1;X(5*P+i)*X(5*P+1:6*P)];                    % thD thD
end
C1 = []; C2 = []; C3 = []; C4 = [];
for i=1:P
    B1 = [B1;X(i)*A3;];               % x Vx Vx
    B2 = [B2;X(i)*A4;];               % x Vy Vy
    B3 = [B3;X(P+i)*A3;];             % y Vx Vx
    B4 = [B4;X(P+i)*A4;];             % y Vy Vy
    D2 = [D2;X(i)*D1];                  % x thD thD
    D3 = [D3;X(i+P)*D1];                % y thD thD
    C1 = [C1,E1(:,(i-1)*P+1:i*P)*X(1:P)];             % X
    C2 = [C2,E1(:,(i-1)*P+1:i*P)*X(P+1:2*P)];         % y
%     C3=[C3,E1_L2(:,(i-1)*P+1:i*P)*X(1:P)];          % X L^2
%     C4=[C4,E1_L2(:,(i-1)*P+1:i*P)*X(P+1:2*P)];      % Y L^2
    %C3=[C3,-E1(:,(i-1)*P+1:i*P)*2*X(1:P)];
    %C4=[C4, E1(:,(i-1)*P+1:i*P)*2*X(P+1:2*P)];
    
end
M = zeros(3*P,3*P);
M(1:P,1:P) = m*E; M(P+1:2*P,P+1:2*P) = m*E; M(2*P+1:3*P,2*P+1:3*P) = m*E_L2;
%M(1:P,1:P)=m*E_L2; M(P+1:2*P,P+1:2*P)=m*E_L2; M(2*P+1:3*P,2*P+1:3*P)=m*E_L4;
Cq1 = [E zeros(P,P);zeros(P,P) E; C2 -C1];
Cq2 = [E zeros(P,P) C2; zeros(P,P) E -C1];
%Cq2=[E_L2 zeros(P,P) C4; zeros(P,P) E_L2 -C3];

term3 = -m*g*E3;
%term3= -m*g*E4;
%term4= -E2*B1 - E2*B2;
term4 = -E2*D2;
%term5= -E2*B3 - E2*B4;
term5 = -E2*D3;

LHS = [M Cq1;Cq2 zeros(2*P,2*P)];
RHS = [zeros(P,1);term3;zeros(P,1);term4;term5];

% CCC = cond(LHS);
% fprintf("time,%d\tcond:%d\n",t,CCC);
Y = linsolve(LHS,RHS);
dxdt(DOF*P+1:2*DOF*P) = Y(1:DOF*P);

