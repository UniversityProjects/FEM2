function z=uexact(x,y,k)
% -------
% a due componenti, la componente (1 oppure 2) e' decisa da k
% -------

if k==1        % prima componente
    z = 0.25*((x^2 - 1)^2)*(y^2 - 1)*y;
elseif k==2    % seconda componente
    z = 0.25*((y^2 - 1)^2)*(1 - x^2)*x;
end
