%
function [uh, errL2, errH1] = main_P2_Vect(xv,yv,vertices,edges,endpoints,boundary,boundedges);
%
% ---------------------------------------------------------------
% FEM per k=2
% Adattamento del codice per la Diffusione di grado 2 al caso
% ... vettoriale. Passaggio preliminare allo Stokes P2/P0.
% Ordinamento globale: prima tutti i GdL della prima componente,
% ... poi tutti i GdL della seconda componente.
% Questa versione fa le matrici locali gia vettoriali e poi
% ... assembla al globale (ordinamento locale: prima GdL della
% ... prima componente, poi della seconda)
% ---------------------------------------------------------------
% addpath ../../Meshes/
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
%
% gradienti delle funzioni di base calcolati
% nei nodi di quadratura dell'elemento di riferimento
%
gphihqx = zeros(6,Nq);
gphihqy = zeros(6,Nq);
%
for i=1:6
    for q=1:Nq
        [gx,gy] = gradhphih2(i,xhq(q),yhq(q));
        gphihqx(i,q) = gx;
        gphihqy(i,q) = gy;
    end
end
%
% ASSEMBLIAMO LA MATRICE DEI COEFFICIENTI
%
A = sparse(2*(nver+nedge),2*(nver+nedge));
b = zeros(2*(nver+nedge),1);
% b1 = zeros(nver+nedge,1);   % parte carico prima componente
% b2 = zeros(nver+nedge,1);   % parte carico seconda componente
%
for iele=1:nele
    %
    % MATRICE ELEMENTARE
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
    KE = zeros(6,6);
    %
    for i=1:6
        for j=1:i-1
            KE(i,j) = KE(j,i);
        end
        %
        for j=i:6
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
                tmp = dot(...
                      (JFIT*[gphihqx(j,q);gphihqy(j,q)]),...
                      (JFIT*[gphihqx(i,q);gphihqy(i,q)])...
                      );
                %
                KE(i,j) = KE(i,j) + c(xq,yq)*tmp*whq(q);
            end
            %
            KE(i,j) = 2*area*KE(i,j);
        end
    end
    KE = [KE,zeros(6);zeros(6),KE]; % matrice duplicata ("vettoriale")
    %
    % ASSEMBLIAMO LA MATRICE GRANDE
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
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];   % prima componente
    dofg = [dofg , dofg + (nver+nedge)];     % aggiungo seconda componente
    %
    A(dofg,dofg) = A(dofg,dofg) + KE;
    %
end

% ---------------------------------
% CALCOLO IL TERMINE DI CARICO 
% Le due componenti fatte separatamente usando P2load.m
% ... e poi combinate. Nota: e' gia ritagliato dei dofs di bordo.
% ---------------------------------

b1 = P2load1(xv,yv,vertices,edges,boundary,boundedges,endpoints);
b2 = P2load2(xv,yv,vertices,edges,boundary,boundedges,endpoints);
b = [b1;b2];

% ---------------------------------
% sistemiamo le condizioni al bordo
% ---------------------------------

internV = setdiff([1:1:nver],boundary); % indici dei vertici INTERNI
internE = setdiff([1:1:nedge],boundedges); % indici dei lati INTERNI
NL = [internV,nver+internE];     % gradi di libert? interni

% Ora considero che i GDL sono per prima e poi seconda componente:
NL2 = (nver+nedge) + NL;  % GdL interni associati alla 2nda componente
NL  = [NL,NL2];           % tutti i GdL interni
clear NL2;

uh = zeros(2*(nver+nedge),1);
%
% estraiamo la sottomatrice corrispondente
% ai nodi liberi NL
%
Kh = A(NL,NL);
fh = b;
%
uh(NL) = Kh\fh;
uh1 = uh(1:length(uh)/2);     % estraggo la prima componente
uh2 = uh(length(uh)/2+1:end);  % estraggo la seconda componente
%
% se iv e' un vertice --> uh(iv) = valore di uh in 
%                                  quel vertice
% se ie e' un edge --> uh(nver+ie) = valore di uh
%                                    nel punto medio
%                                    dell'edge
%
% Adesso calcoliamo l'errore in L2 rispetto a ue
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
errL2sq = 0;   % errore L2
errH1sq = 0;   % errore H1
uH1 = 0;       % norma H1 della sol. esatta (serve per err. relativo)
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
      sq = sq + (uexact(xq,yq,1)-tmp1)^2*whq(q) + (uexact(xq,yq,2)-tmp2)^2*whq(q);     
      sqH = sqH + norm(uexactG(xq,yq,1) - tmpH1)^2*whq(q) ...
                + norm(uexactG(xq,yq,2) - tmpH2)^2*whq(q);
      uH = uH + (norm(uexactG(xq,yq,1))^2+norm(uexactG(xq,yq,2))^2)*whq(q);       
    end
    %
    sq = sq*2*area;
    sqH = sqH*2*area;
    uH = uH*2*area;
    %
    errL2sq = errL2sq + sq;
    errH1sq = errH1sq + sqH;
    uH1 = uH1 + uH;
    %
end
%
errL2 = sqrt(errL2sq)  % (comprende ambo le componenti)
errH1 = sqrt(errH1sq/uH1)  % (comprende ambo le componenti) 
%


% % ----------------------------------
% % PLOT DELLA SOLUZIONE uh COMPONENTE PER COMPONENTE
% % (basato sui valori ai vertici)
% % ----------------------------------
% figure(1)
% for k=1:size(vertices,1)
%     hold on;
%     index=(vertices(k,1:3))';
%     % riga sotto: opzione per funzione nota ai vertici
%     uu_tmp=[uh1(index);uh1(index(1))];  % prima componente
%     % uu_tmp=[uh2(index);uh2(index(1))];  % seconda componente
%     % riga sotto: opzione per funzione costante a tratti
%     % p_tmp = ph(k)*ones(length(index)+1,1); 
%     vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
%     fill3(vert_temp(:,1),vert_temp(:,2),uu_tmp,uu_tmp); 
% end
% view(3)
% grid on
% colorbar
% hold off        
        
% ----------------------------------
% PLOT DELLA SOLUZIONE uh COME FRECCE AI VERTICI
% ----------------------------------        
%       
% figure(2)
% quiver(xv,yv,uh1(1:nver),uh2(1:nver));
%
%