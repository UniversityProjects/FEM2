% ----------------------------------------------------
function []=meshplot(xv,yv,endpoints);
% -----------------------------------------------------------------
% given some mesh data (see meshgen for format), it plots the mesh.
% -----------------------------------------------------------------
figure(1)
hold on
P = [xv,yv];
E = endpoints;
for i=1:size(E,1)
    x=[P(E(i,1),1) , P(E(i,2),1)];
    y=[P(E(i,1),2) , P(E(i,2),2)];
    % if E(i,4)==0
    %   plot(x,y,'r')
    %else
       plot(x,y)
    %end    
end
