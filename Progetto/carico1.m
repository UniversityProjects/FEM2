function f1 = carico1(x,y)
%CARICO1
%    F1 = CARICO1(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    17-Feb-2017 22:16:57

t3 = x.^2;
t2 = t3-1.0;
t4 = y.^2;
t5 = t4-1.0;
f1 = t2.^2.*y.*(-3.0./2.0)-t2.*t5.*y-t3.*t5.*y.*2.0;
