function [u]= exsol(x,y);
% date le coordinate globali X=(x,y) 
% ...calcola la soluzione esatta 

% -----------------
% PROBLEMA STAZIONARIO (alfa=1,gamma=0)
% dominio quadrato [0,1]^2
% u = sin(x*pi).*sin(y*pi); 
% -----------------

% -----------------
% PROBLEMA STAZIONARIO (alfa=1,gamma=0)
% su dominio circolare centrato nell'origine, raggio 1
u = 1 - (x.^2 + y.^2); 
% -----------------