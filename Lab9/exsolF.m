function [fex]=exsolF(X);
% date le coordinate globali X=(x,y) calcola la sol esatta (flussi)
% ... in quel punto
x=X(1);
y=X(2);
fex = zeros(2,1);

% -----------------
% TEST 1 su QUADRATO UNITARIO 
fex(1) = (1-2*x)*y*(1-y);
fex(2) = x*(1-x)*(1-2*y);
