function [u]= exsolG(x,y);
% date le coordinate globali X=(x,y) 
% ...calcola il gradiente della soluzione esatta al tempo finale T
% -----------------


% -----------------
% PROBLEMA PARABOLICO UNO
%
% (Il valore gamma e' un fattore che rende la soluzione più rapidamente
% ... variabile in tempo)
%
T = 1;
gamma = 1;
u = exp(T*gamma)*pi*[cos(pi*x).*sin(pi*y);sin(pi*x).*cos(pi*y)];
% u = exp(T*gamma)*(sin(pi*x).*sin(pi*y));  <--- soluzione
% -----------------



