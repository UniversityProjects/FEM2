function z=uexact(x,y,k)
% -------
% a due componenti, la componente (1 oppure 2) e' decisa da k
% -------

% PROBLEMA "TRIGONOMETRICO" (STOKES)

if k==1        % prima componente
    z = -(cos(pi*(x-0.5))^2)*cos(pi*(y-0.5))*sin(pi*(y-0.5))/2;
elseif k==2    % seconda componente
    z = (cos(pi*(y-0.5))^2)*cos(pi*(x-0.5))*sin(pi*(x-0.5))/2;
end

