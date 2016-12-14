function [u]= exsolG(x,y);
% date le coordinate globali X=(x,y) 
% ...calcola il GRADIENTE della soluzione esatta 

% -----------------
% PROBLEMA STAZIONARIO (alfa=1,gamma=0)
% domninio quadrato [0,1]^2
% u = pi*[cos(x*pi).*sin(y*pi);sin(x*pi).*cos(y*pi)];
% -----------------

% -----------------
% PROBLEMA STAZIONARIO (alfa=1,gamma=0)
% su dominio circolare centrato nell'origine, raggio 1
u = (-2)*[x;y]; 
% -----------------