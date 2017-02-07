function [] = carico (lambda) 

syms x y real;


% Prima Componente
ue = 0.25*((x^2 - 1)^2)*(y^2 - 1)*y;

uex = diff(ue,x);
uey = diff(ue,y);
f = -diff(uex,x) -diff(uey,y) - lambda*(uex + uey);
car1 = matlabFunction(f,'File','carico1');


% Seconda Componente
ue = 0.25*((y^2 - 1)^2)*(1 - x^2)*x;

uex = diff(ue,x);
uey = diff(ue,y);
f = -diff(uex,x) -diff(uey,y) - lambda*(uex + uey);
car2 = matlabFunction(f,'File','carico2');

reset(symengine)