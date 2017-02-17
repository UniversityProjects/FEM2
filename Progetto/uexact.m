function z=uexact(x,y,k)
% -------
% a due componenti, la componente (1 oppure 2) e' decisa da k
% -------

% PROBLEMA "TRIGONOMETRICO"
%if k==1        % prima componente
%    z = sin(x*pi).*sin(y*pi);
%elseif k==2    % seconda componente
%    z = 3*sin(x*pi).*sin(y*pi);
%end

if k==1        % prima componente
    z = 0.25*((x^2 - 1)^2)*(y^2 - 1)*y;
    % z = cos((pi/2)*x)*cos((pi/2)*y);
elseif k==2    % seconda componente
    z = 0.25*((y^2 - 1)^2)*(1 - x^2)*x;
    % z = cos((pi/2)*x)*cos((pi/2)*y);
end
