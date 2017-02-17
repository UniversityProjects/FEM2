function [] = carico (lambda) 

syms x y real;

% Prima Componente
ue1 = 0.25*((x^2 - 1)^2)*(y^2 - 1)*y;
% ue1 = cos((pi/2)*x)*cos((pi/2)*y);

% Seconda Componente
ue2 = 0.25*((y^2 - 1)^2)*(1 - x^2)*x;
% ue2 = cos((pi/2)*x)*cos((pi/2)*y);

% Derivate Prime Prima Componente
ue1x = diff(ue1,x);
ue1y = diff(ue1,y);

% Derivate Prime Seconda Componente
ue2x = diff(ue2,x);
ue2y = diff(ue2,y);

% Calcolo Termine Noto
f1 = -diff(ue1x,x) - diff(ue1y,y) - lambda*diff((ue1x + ue2y),x);
f2 = -diff(ue2x,x) - diff(ue2y,y) - lambda*diff((ue1x + ue2y),y);

% Scrittura Su File Delle Componenti Del Termine Noto
car1 = matlabFunction(f1,'File','carico1');
car2 = matlabFunction(f2,'File','carico2');


% Seconda Componente






reset(symengine)