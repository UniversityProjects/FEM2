function [f]= carico(x,y,t);
% date le coordinate globali X=(x,y) e il tempo t
% ...calcola il carico (scalare) in quel punto e instante temporale


% -----------------
% PROBLEMA PARABOLICO UNO
%
% (Il valore gamma e' un fattore che rende la soluzione più rapidamente
% ... variabile in tempo)
%
gamma = 1;
f = exp(gamma*t)*(2*(pi^2)+gamma)*(sin(pi*x)*sin(pi*y));
% -----------------


