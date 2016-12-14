function [gx,gy]=gradhphih2(i,x,y)
%
% gradienti delle funzioni di base sull'elemento di riferimento
% per l'elemento MINI. 
%
switch i
    case 1
        gx = -1;
        gy = -1;
    case 2
        gx = 1;
        gy = 0;
    case 3
        gx = 0;
        gy = 1;
    case 4                 % bolla interna (vale 1 nel baricentro)
        gx = 27*(y-2*x*y-y^2);
        gy = 27*(x-x^2-2*x*y);
end
