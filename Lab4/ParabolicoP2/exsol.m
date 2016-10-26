function [u]= exsol(x,y);
% date le coordinate globali X=(x,y) 
% ...calcola la soluzione esatta al tempo finale T

% -----------------
% PROBLEMA PARABOLICO UNO
%
% (Il valore gamma e' un fattore che rende la soluzione più rapidamente
% ... variabile in tempo)
%
T = 1;
gamma = 10;
u = exp(gamma*T)*(sin(pi*x).*sin(pi*y));
% -----------------



