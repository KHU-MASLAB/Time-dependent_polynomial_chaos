function [Condition, isterminal, direction] = stopopt(te,y,t)
global stopcrt_t
%Condition = abs(max(y(2:P))) < abs(y(1))/stopcrt;
Condition = te > stopcrt_t;
isterminal = 1; 
direction = 0;
