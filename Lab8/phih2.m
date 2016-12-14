function z = phih2(i,xh,yh)
% 
% funzioni di base sull'elemento di riferimento
% per l'elemento MINI. 
%
switch i
    case 1
        z = 1 - xh - yh;
    case 2
        z = xh;
    case 3
        z = yh;
    case 4                 % bolla interna (vale 1 nel baricentro)
        z = 27*xh*yh*(1-xh-yh);
end
