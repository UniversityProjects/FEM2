function [f]= carico(X);
% date le coordinate globali X=(x,y) calcola il carico (scalare)
% ... in quel punto
x=X(1);
y=X(2);

% -----------------
% TEST 1 su QUADRATO UNITARIO 
f = 2*y*(1-y) + 2*x*(1-x);
% -----------------
