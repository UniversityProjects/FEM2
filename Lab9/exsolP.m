function [pex]=exsolP(X);
% date le coordinate globali X=(x,y) calcola la sol esatta (pressioni)
% ... in quel punto
x=X(1);
y=X(2);

% -----------------
% TEST 1 su QUADRATO UNITARIO 
pex = x*(1-x)*y*(1-y);