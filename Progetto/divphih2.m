function [z]=divphih2(i,x,y)

% Divergenza come somma delle componenti del gradiente
[gx, gy] = gradhphih2 (i, x, y);

z = gx + gy;