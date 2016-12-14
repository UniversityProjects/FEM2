function z=carico(x,y,k)
% -------
% carico a due componenti, la componente (1 oppure 2) e' decisa da k
% -------

% PROBLEMA "TRIGONOMETRICO" (STOKES)

if k==1        % prima componente
 z = (pi^2)*sin(pi*(y-0.5))*cos(pi*(y-0.5))*(1-4*(cos(pi*(x-0.5)))^2) - pi*cos(pi*(x-0.5));
elseif k==2    % seconda componente
 z = (-pi^2)*sin(pi*(x-0.5))*cos(pi*(x-0.5))*(1-4*(cos(pi*(y-0.5)))^2) + pi*cos(pi*(y-0.5));
end

