function [XDot,acc] = SolveDAE_stochastic(t,X,Prop,Galerkin)

P = length(Galerkin.E2);
RDof = 3;

% rigid body
q1 = X(1:RDof*P);
q1Dot = X(RDof*P+1:(2*RDof)*P);

%% Constraint
[phiq, gamma] = Constraint_PC(q1,q1Dot,Galerkin); 

%% mass matirx
M = Mass_PC(Prop,Galerkin);

%% Force Vector
Qext = kron([0; -Prop.m1*Prop.g; 0],Galerkin.E1);

%% Assemble
dimLambda = size(phiq,1)/P;
MM = [M phiq' ; phiq  zeros(dimLambda*P,dimLambda*P)];
Q = [Qext ; gamma];

%% Solve System
acc = linsolve(MM,Q);

qDotDot = acc(1:RDof*P,:);
lambda = acc(RDof*P+1:(RDof+dimLambda)*P,:); 

XDot = [q1Dot ; qDotDot];


end







