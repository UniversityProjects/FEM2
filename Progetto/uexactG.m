function Z=uexactG(x,y,k)
% -------
% considera due componenti, la componente (1 oppure 2) e' decisa da k
% calcola gradiente della sol esatta
% -------

% PROBLEMA "TRIGONOMETRICO"

% if k==1        % prima componente
%     Z = pi*[cos(x*pi).*sin(y*pi);sin(x*pi).*cos(y*pi)];
% elseif k==2    % seconda componente
%     Z = 3*pi*[cos(x*pi).*sin(y*pi);sin(x*pi).*cos(y*pi)];
% end


if k==1        % prima componente
    %Z = x*y*(x^2 - 1)*(y^2 - 1);
    Z = -(pi/2)*sin((pi/2)*x)*cos((pi/2)*y);
elseif k==2    % seconda componente
    %Z = -x*y*(x^2 - 1)*(y^2 - 1);
    Z = -(pi/2)*cos((pi/2)*x)*sin((pi/2)*y);
end

