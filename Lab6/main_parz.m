%
function [uh,ph] = main_parz(xv,yv,vertices,edges,endpoints,boundary,boundedges);
%
% ---------------------------------------------------------------
% FEM per k=2
% Risolve Stokes con il P2P0. Usa il grad-grad (non simmetrico).
%
% Ordinamento globale: prima tutti i GdL della prima componente,
% ... (vertici,lati) e poi tutti i GdL della 2a componente.
% Ordinamento locale: prima componente nei tre vertici, poi 
% ... prima componente nei tre lati, poi 2a componente nei tre
% ... vertici, poi 2a nei tre lati
% NOTA: la parte A della matrice globale e' fatta per duplicazione
% ... diretta dal caso scalare (anche al locale),
% ... mentre la parte B e' costruita gia correttamente al locale
% ---------------------------------------------------------------
% NOTA: manca il calcolo dell'errore per la pressione
% ---------------------------------------------------------------
%

% clear all
% close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh Load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[xv,yv,vertices,edges,endpoints,boundary,boundedges] = load('meshesquad.mat');
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
B = sparse(2*(nver+nedge),nele);
b1 = zeros(nver+nedge,1);   % parte carico prima componente
b2 = zeros(nver+nedge,1);   % parte carico seconda componente
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
    VettArea(iele)=area;  % questo vettore delle aree servira' in seguito,
                          % ... e' un vettore che ha in posizione 
                          % ... i-esima l'area del i-esimo elemento
    % ------------
    KE = zeros(6,6);    % matrici locali per la A (vedi NOTA in cima)
    BE = zeros(12,1);   % matrici locali per la B (in realta' un vettore)
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
    
    % ----------------------
    % parte per la DIVERGENZA (matrice B)
    % (si calcola l'integrale della DIV sul elemento, poiche' pressioni P0)
    % ----------------------
    %        
    for j=1:12  % 12 e' il numero di GdL locali considerando che e' vettoriale
            for q=1:Nq
                %
                %
                % tmp = JF*[xhq(q);yhq(q)]+[x1;y1];
                %
                % xq = tmp(1);
                % yq = tmp(2);
                %
      % NOTA: la DIV delle prime sei funzioni di base (associate alla 1a componente)
      % ... e' data dalla loro derivata rispetto ad x, cioe' la prima componente
      % del loro vettore gradiente. 
      % Analogamente per le ultime sei (2a componente) e' la seconda
      % ... componente del loro gradiente.
      % Si deve anche trasformare il gradiente sul elemento fisico
      % ... prima di selezionare la componente.
            if j < 6.5         % prime sei funzioni di base
                tmp = JFIT(1,1:2)*[gphihqx(j,q);gphihqy(j,q)];
            elseif j > 6.5   % seconde sei funzioni di base
                tmp = JFIT(2,1:2)*[gphihqx(j-6,q);gphihqy(j-6,q)];
            end
            %
            BE(j)= BE(j) + tmp*whq(q);                             
            end
       BE(j) = 2*area*BE(j);     
    end
    
    % ------------------------------------------------------
    % ASSEMBLIAMO LA MATRICE A (versione scalare, poi va "duplicata")
    % ------------------------------------------------------
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
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];   % 1a componente
    dofgg = [dofg , dofg + (nver+nedge)];    % aggiungo anche 2a componente    
    %   
    A(dofgg,dofgg) = A(dofgg,dofgg) + KE;     % assemblo
    %
    % termine noto
    %
    fE1 = zeros(6,1);   % carico per la prima componente
    fE2 = zeros(6,1);   % carico per la seconda componente
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
            fE1(i) = fE1(i) + carico(xq,yq,1)*phihq(i,q)*whq(q);
            fE2(i) = fE2(i) + carico(xq,yq,2)*phihq(i,q)*whq(q);
        end
        %
        fE1(i) = 2*area*fE1(i);
        fE2(i) = 2*area*fE2(i);
        %
    end
    %
    b1(dofg) = b1(dofg) + fE1;
    b2(dofg) = b2(dofg) + fE2;
    %
   % ------------------------------------------------------
   % ASSEMBLIAMO LA MATRICE B
   % ------------------------------------------------------
   B(dofgg,iele) = B(dofgg,iele) + BE;     
end

b = [b1;b2];    % combino le due componenti del carico in singolo vettore

% ---------------------------------
% sistemiamo le condizioni al bordo
% ---------------------------------

internV = setdiff([1:1:nver],boundary); % indici dei vertici INTERNI
internE = setdiff([1:1:nedge],boundedges); % indici dei lati INTERNI
NL = [internV,nver+internE];     % gradi di liberta' interni

% Ora considero che i GDL sono per prima e poi seconda componente:
NL2 = (nver+nedge) + NL;  % GdL interni associati alla 2nda componente
NL  = [NL,NL2];           % tutti i GdL interni
clear NL2;

%
% estraiamo la sottomatrice corrispondente
% ai nodi liberi NL
%

Ah = A(NL,NL);
Bh = B(NL,:);
fh = b(NL);
%
clear A
clear B
clear b
% -----------------------------------------------------
% Costruisco la matrice totale (parti A e B) e il carico
%
Kh = [Ah, Bh; Bh', zeros(nele)];
fh = [fh; zeros(nele,1) ];

% -----------------------------------------------------
% SE le condizioni al bordo sulla velocita' sono tutte
% ... Dirichlet (omogenee), si deve imporre condizione 
% ... sulla pressione (media nulla). Lo facciamo "con moltiplicatore".
% ALTRIMENTI commentare le TRE righe sotto.
% -----------------------------------------------------
N = length(NL);
Kh = [Kh, [zeros(N,1); VettArea']; [zeros(1,N), VettArea, 0] ];
fh = [fh; 0];


uh = zeros(2*(nver+nedge),1);
% -------------------------------
% Risolviamo il sistema lineare
% -------------------------------
solh = Kh\fh;
uh(NL) = solh(1:N);            % estraggo le velocità
uh1 = uh(1:length(uh)/2);      % estraggo la prima componente
uh2 = uh(length(uh)/2+1:end);  % estraggo la seconda componente
ph  = solh(N+1:N+nele);        % estraggo le pressioni


% -------------------------------------------------------
% Adesso calcoliamo l'errore in L2 e H1 rispetto a ue,pe
% -------------------------------------------------------
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
errL2sq = 0;   % errore L2 (velocità)
errH1sq = 0;   % errore H1 (velocità)
preL2sq = 0;   % errore L2 (pressioni) 
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
    end
    %
    sq = sq*2*area;
    sqH = sqH*2*area;
    %
    errL2sq = errL2sq + sq;
    errH1sq = errH1sq + sqH;
    %
end
%
disp('Errore in L2 delle velocità')
errL2 = sqrt(errL2sq)  % (comprende ambo le componenti)
disp('Errore in H1 delle velocità')
errH1 = sqrt(errH1sq)  % (comprende ambo le componenti) 
%


% ----------------------------------
% PLOT DELLA SOLUZIONE uh COMPONENTE PER COMPONENTE
% (basato sui valori ai vertici)
% ----------------------------------
figure(1)
for k=1:size(vertices,1)
    hold on;
    index=(vertices(k,1:3))';
    % riga sotto: opzione per funzione nota ai vertici
    uu_tmp=[uh1(index);uh1(index(1))];  % prima componente
    % uu_tmp=[uh2(index);uh2(index(1))];  % seconda componente
    vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
    fill3(vert_temp(:,1),vert_temp(:,2),uu_tmp,uu_tmp); 
end
view(3)
grid on
colorbar
hold off   

% ----------------------------------
% PLOT DELLA SOLUZIONE ph 
% ----------------------------------

figure(2)
for k=1:nele
    hold on;
    index=(vertices(k,1:3))';
    % riga sotto: opzione per funzione costante a tratti
    p_tmp = ph(k)*ones(length(index)+1,1); 
    vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
    fill3(vert_temp(:,1),vert_temp(:,2),p_tmp,p_tmp); 
end
view(3)
grid on
colorbar
hold off   

% ----------------------------------
% PLOT DELLA SOLUZIONE uh COME FRECCE AI VERTICI
% ----------------------------------        
%       
figure(3)
quiver(xv,yv,uh1(1:nver),uh2(1:nver));
%
%