function z=uexact(x,y,k)
% -------
% a due componenti, la componente (1 oppure 2) e' decisa da k
% -------

% PROBLEMA "TRIGONOMETRICO"

if k==1        % prima componente
    z = sin(x*pi).*sin(y*pi);
elseif k==2    % seconda componente
    z = 3*sin(x*pi).*sin(y*pi);
end

