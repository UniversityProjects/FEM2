function Z=uexactG(x,y,k)
% -------
% considera due componenti, la componente (1 oppure 2) e' decisa da k
% calcola gradiente della sol esatta
% -------

% ------------------------------------------
% PROBLEMA "TRIGONOMETRICO"
%
x = pi*(x-0.5);
y = pi*(y-0.5);
if k==1        % prima componente
    Z = pi*[cos(x)*sin(x)*cos(y)*sin(y); -1/2*(cos(x)^2)*(-sin(y)^2 + cos(y)^2)];
elseif k==2    % seconda componente
    Z = pi*[1/2*(cos(y)^2)*(-sin(x)^2 + cos(x)^2);-cos(x)*sin(x)*cos(y)*sin(y)];
end
% -----------------------------------------


