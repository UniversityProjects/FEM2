function b = P2load(xv,yv,vertices,edges,boundary,boundedges,endpoints);
%
% Modifica del file dei P2 che hanno gia da Ale Russo.
% Calcola il vettore termine noto (dato tramite carico.m) al tempo t.
%

%
% FORMULA DI QUADRATURA 
%
fdq = 'degree=2';
%
[xhq,yhq,whq]=quadratura(fdq);
%
Nq = length(xhq);
%
nver  = length(xv);         % calcolo numero vertici
nedge = size(endpoints,1);  % calcolo numero lati
nele = size(vertices,1);    % calcolo numero elementi
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
    end
end

% ----------------

b = zeros(nver+nedge,1); 

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
    
    % termine noto
    %
    fE = zeros(6,1);
    %
    for i=1:6
        for q=1:Nq
            %
            % calcoliamo l'immagine (xq,yq) sul triangolo
            % corrente del nodo di quadratura
            % (xhq(q),yhq(q)) che sta
            % sull'elemento di rif.
            %
            tmp = JF*[xhq(q);yhq(q)]+[x1;y1];
            %
            xq = tmp(1);
            yq = tmp(2);
            %
            fE(i) = fE(i) + carico(xq,yq)*phihq(i,q)*whq(q);
        end
        %
        fE(i) = 2*area*fE(i);
        %
    end
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
    % array dei dof globali
    %
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];
    %
    b(dofg) = b(dofg) + fE;
    %
end

% --------------------------------------------------
% Adesso ritagliamo il carico ai gradi di libertà interni 
% (cond al bordo: Dirichlet omogenee)
% --------------------------------------------------

nver = length(xv);          % calcolo numero vertici
nedge = size(endpoints,1);  % calcolo numero lati
internV = setdiff([1:1:nver],boundary); % indici dei vertici INTERNI
internE = setdiff([1:1:nedge],boundedges); % indici dei lati INTERNI
NL = [internV,nver+internE];     % gradi di libertà interni

b = b(NL);



    
    
    
    











        
        
        




















    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





























