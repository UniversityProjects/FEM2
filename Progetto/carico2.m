function f2 = carico2(x,y)
%CARICO2
%    F2 = CARICO2(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    23-Feb-2017 22:49:05

t3 = y.^2;
t2 = t3-1.0;
t4 = x.^2;
t5 = t4-1.0;
f2 = t2.^2.*x.*(3.0./2.0)+t2.*t5.*x+t3.*t5.*x.*2.0;
