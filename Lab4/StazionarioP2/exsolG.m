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
% u = (-2)*[x;y]; 
% -----------------

% -----------------
% PROBLEMA STAZIONARIO (alfa=1,gamma=0)
% su dominio (poligonale) a forma di boomerang, 
% vedi dentro cartella "Meshes/saved_meshes"
% Caricare il file boomerang_mesh.mat
%
u = [-((x - 1).*(- 2*x.^2 + 4*x + 5*y.^2 - 6*y))/2 ; ...
    -(5*(x.^2).*y)/2 + (3*x.^2)/2 + 5*x.*y - 3*x + 4*y.^3 - 9*y.^2 + 4*y];

% -----------------