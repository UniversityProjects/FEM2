function z=pexact(x,y)
%
% pressione esatta (soluzione)
%
% NOTA: nulla nel vertice di coord (0,0), ed e' anche nulla la media.

% PROBLEMA "TRIGONOMETRICO" (STOKES)

z = (sin(pi*(x-0.5)) - sin(pi*(y-0.5)));
