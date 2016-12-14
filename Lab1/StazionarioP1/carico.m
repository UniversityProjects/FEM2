function [f]= carico(x,y);
% date le coordinate globali X=(x,y) e il tempo t
% ...calcola il carico (scalare) in quel punto e instante temporale


% -----------------
% PROBLEMA PARABOLICO UNO
% f = exp(t)*(2*(pi^2)+1)*(sin(pi*x)*sin(pi*y));
% -----------------

% -----------------
% PROBLEMA STAZIONARIO (alfa=1,gamma=0)
% dominio quadrato [0,1]^2
% f = sin(x*pi).*sin(y*pi)*(2*pi^2);
% -----------------

% -----------------
% PROBLEMA STAZIONARIO (alfa=1,gamma=0)
% su dominio circolare centrato nell'origine, raggio 1
f = 4; 
% -----------------