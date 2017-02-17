function [errL2, errH1] = err(uh, fdq, xv, yv, vertices, edges, endpoints, out)

uh1 = uh(1:length(uh)/2);     % estraggo la prima componente
uh2 = uh(length(uh)/2+1:end);  % estraggo la seconda componente

nver  = length(xv);         % calcolo numero vertici
nedge = size(endpoints,1);  % calcolo numero lati
nele = size(vertices,1);    % calcolo numero elementi

%
% FORMULA DI QUADRATURA PER CALCOLARE L'ERRORE
%
fdq = 'degree=5';
%
[xhq,yhq,whq]=quadratura(fdq);
%
Nq = length(xhq);
%
% (xhq,yhq) = nodi di quadratura su T cappello
% whq = relativi pesi
%
% funzioni di base calcolate nei nodi di 
% quadratura dell'elemento di riferimento
%
phihq = zeros(6,Nq);
%
for i=1:6
    for q=1:Nq
        phihq(i,q) = phih2(i,xhq(q),yhq(q));
        [gx,gy] = gradhphih2(i,xhq(q),yhq(q));
        gphihqx(i,q) = gx;
        gphihqy(i,q) = gy;
    end
end
%
errL2sq = 0;    % L2 error
errH1sq = 0;    % H1 error
uL2 = 0;        % norma L2 della sol. esatta (serve per err. relativo)
uH1 = 0;        % norma H1 della sol. esatta (serve per err. relativo)
%
for iele=1:nele
    %
    v1 = vertices(iele,1);
    v2 = vertices(iele,2);
    v3 = vertices(iele,3);
    %
    x1 = xv(v1);
    y1 = yv(v1);
    %
    x2 = xv(v2);
    y2 = yv(v2);
    %
    x3 = xv(v3);
    y3 = yv(v3);
    %
    % Jacobiana della trasformazione F
    %
    JF = [x2-x1 x3-x1
          y2-y1 y3-y1];
    %
    % inversa
    %
    JFI = inv(JF);
    %
    % trasposta
    %
    JFIT = JFI';
    %
    area = 1/2*det(JF);
    %
    % recuperiamo i coefficienti delle phi_i
    %
    % gradi di liberta' globali:
    %
    % vertice i  ----> i
    %
    % lato l -----> nver+l
    %
    % recuperiamo i lati del triangolo
    %
    l1 = edges(iele,1);
    l2 = edges(iele,2);
    l3 = edges(iele,3);
    %
    % array dei dof globali del triangolo T
    %
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];
    %
    uT1 = uh1(dofg);   % prima componente
    uT2 = uh2(dofg);   % seconda componente 
    %
    sq  = 0;
    sqH = 0;
    uL = 0;
    uH  = 0;
    %
    for q=1:Nq
        %
        % calcolo la sommatoria sulle phi_i
        %
        tmp1  = 0;       % valori della funzione uh (prima componente)
        tmp2  = 0;       % (seconda componente)
        tmpH1 = [0;0];   % valori dei gradienti di uh (prima componente)
        tmpH2 = [0;0];   % (seconda componente)
        %
        for i=1:6
            tmp1 = tmp1 + uT1(i)*phihq(i,q);
            tmp2 = tmp2 + uT2(i)*phihq(i,q);
            tmpH1= tmpH1 + uT1(i)*(JFIT*[gphihqx(i,q);gphihqy(i,q)]); 
            tmpH2= tmpH2 + uT2(i)*(JFIT*[gphihqx(i,q);gphihqy(i,q)]);
        end
        %
        pos = JF*[xhq(q);yhq(q)]+[x1;y1];
        %
        xq = pos(1);
        yq = pos(2);
        %
%       sq = sq + (uexact(xq,yq,1)-tmp1)^2*whq(q) + (uexact(xq,yq,2)-tmp2)^2*whq(q);     
%       sqH = sqH + norm(uexactG(xq,yq,1) - tmpH1)^2*whq(q) ...
%                 + norm(uexactG(xq,yq,2) - tmpH2)^2*whq(q);
%       uL = uL + (norm(uexact(xq,yq,1))^2+norm(uexact(xq,yq,2))^2)*whq(q);
%       uH = uH + (norm(uexactG(xq,yq,1))^2+norm(uexactG(xq,yq,2))^2)*whq(q);
      sq = sq + (uexact(xq,yq,1)-tmp1)^2*whq(q) + (uexact(xq,yq,2)-tmp2)^2*whq(q);     
      sqH = sqH + norm([uexactGmine(xq,yq,1,1); uexactGmine(xq,yq,2,1)] - tmpH1)^2*whq(q) ...
                + norm([uexactGmine(xq,yq,1,2); uexactGmine(xq,yq,2,2)] - tmpH2)^2*whq(q);
      uL = uL + (norm(uexact(xq,yq,1))^2+norm(uexact(xq,yq,2))^2)*whq(q);
      uH = uH + (norm([uexactGmine(xq,yq,1,1); uexactGmine(xq,yq,2,1)])^2 + ...
                 norm([uexactGmine(xq,yq,1,2); uexactGmine(xq,yq,2,2)])^2)*whq(q);       
    end
    %
    sq = sq*2*area;
    sqH = sqH*2*area;
    uL = uL*2*area;
    uH = uH*2*area;
    %
    errL2sq = errL2sq + sq;
    errH1sq = errH1sq + sqH;
    uL2 = uL2 + uL;
    uH1 = uH1 + uH;
    %
end
%

errL2 = sqrt(errL2sq/uL2);  % Relative L2 Error
errH1 = sqrt(errH1sq/uH1);  % Relative H1 Error

if (strcmp(out,'yes'))
    disp(['      L2 Error (Relative): ' num2str(errL2)]);
    disp(['      H1 Error (Relative): ' num2str(errH1)]);
end