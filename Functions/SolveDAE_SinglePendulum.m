function [XDot,acc] = SolveDAE_SinglePendulum(t,X,Prop)

RDof = 3;

% rigid body
q1 = X(1:RDof);
q1Dot = X(RDof+1:2*RDof);

%% Constraint
[phiq, gamma] = Constraint(q1,q1Dot,Prop.L); 

%% mass matirx
M = Prop.mass;

%% Force Vector
Qext = [0; -Prop.m1*Prop.g; 0];

%% Assemble
dimLambda = size(phiq,1);
MM = [M phiq' ; phiq  zeros(dimLambda,dimLambda)];
Q = [Qext ; gamma];

%% Solve System
acc = linsolve(MM,Q);

qDotDot = acc(1:RDof,:);
lambda = acc(RDof+1:RDof+dimLambda,:); 

XDot = [q1Dot ; qDotDot];


end







