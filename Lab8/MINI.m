%
function [errPR, errH1, uh,ph] = MINI(xv,yv,vertices,edges,endpoints,boundary,boundedges);
%
% ---------------------------------------------------------------
% FEM per k=2
% Risolve Stokes con il MINI. Usa il grad-grad (non simmetrico).
%
% Ordinamento globale: prima tutti i GdL della prima componente,
% ... (vertici,lati) e poi tutti i GdL della 2a componente.
% Ordinamento locale: prima componente nei tre vertici, poi 
% ... prima componente della bolla, poi 2a componente nei tre
% ... vertici, poi 2a della bolla
% NOTA: la parte A della matrice globale e' fatta per duplicazione
% ... diretta dal caso scalare (anche al locale),
% ... mentre la parte B e' costruita gia correttamente al locale
% PRESSIONI: i gradi di lib sono associati ai vertici, dunque ordinati come
% ... sono ordinati i vertici (e' un campo scalare).
% ---------------------------------------------------------------
% FUNZIONA SOLO CON COND AL BORDO SU VELOCITA' NULLE, IN QUANTO LA
% MATRICE B VIENE CALCOLATA ANDANDO PER PARTI. SE LO SI VUOLE PIU
% GENERALE "BASTA" MODIFICARE IL CALCOLO DELLA B
% ---------------------------------------------------------------
%
% FORMULA DI QUADRATURA 
%
fdq = 'degree=4';   % piu alta per integrare grad-bolla vs grad-bolla (deg.4)
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
phihq = zeros(4,Nq);
%
for i=1:4
    for q=1:Nq
        phihq(i,q) = phih2(i,xhq(q),yhq(q));
    end
end
%
% gradienti delle funzioni di base calcolati
% nei nodi di quadratura dell'elemento di riferimento
%
gphihqx = zeros(4,Nq);
gphihqy = zeros(4,Nq);
%
for i=1:4
    for q=1:Nq
        [gx,gy] = gradhphih2(i,xhq(q),yhq(q));
        gphihqx(i,q) = gx;
        gphihqy(i,q) = gy;
    end
end
%
% ASSEMBLIAMO LA MATRICE DEI COEFFICIENTI
%
A = sparse(2*(nver+nele),2*(nver+nele));
B = sparse(2*(nver+nele),nver);
b1 = zeros(nver+nele,1);   % parte carico prima componente
b2 = zeros(nver+nele,1);   % parte carico seconda componente
%
vettarea = zeros(nver,1);  % serve per imporre media nulla su ph (vedi sotto)
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
    vettarea([v1 v2 v3]) = vettarea([v1 v2 v3]) + (area/3)*ones(3,1);
    % vettarea e' un vettore che ad ogni vertice associa l'area di 
    % ... tutti gli elementi attorno diviso 3, la quale a sua volta
    % ... corrisponde al integrale della associata funzione di base
    % ... dei P1.
    %
    KE = zeros(4,4);    % matrici locali per la A (vedi NOTA in cima)
    BE = zeros(8,3);    % matrici locali per la B 
    %
    for i=1:4
        for j=1:i-1
            KE(i,j) = KE(j,i);
        end
        %
        for j=i:4
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
    
    KE = [KE,zeros(4);zeros(4),KE]; % matrice duplicata ("vettoriale")
    
    % ----------------------
    % parte per la DIVERGENZA (matrice B)
    % (ricordando che le pressioni sono continue, la calcoliamo come
    % ... gradiente-delle-pressioni prodotto scalare la velocit?)
    % ----------------------
    %        
    for i=1:8  % 8 e' il numero di GdL locali della velocit? (e' vettoriale)
    for j=1:3  % 3 e' il numero di GdL della pressione
            for q=1:Nq
                % uso per semplicita' un ciclo IF per distingure le prime
                % ... quattro (i=1,..,4) che sono associate alla prima 
                % ... componente della velocita', e le ultime quattro 
                % ... (i=5,..,8) che sono associate alla seconda
            if i<4.5 
                tmp = -(JFIT(1,1:2)*[gphihqx(j,q);gphihqy(j,q)])*...
                    phihq(i,q);             
            elseif i > 4.5
                tmp = -(JFIT(2,1:2)*[gphihqx(j,q);gphihqy(j,q)])*...
                    phihq(i-4,q);
            end
            %
            BE(i,j)= BE(i,j) + tmp*whq(q);                             
            end
       BE(i,j) = 2*area*BE(i,j);     
    end
    end
    
    % ------------------------------------------------------
    % ASSEMBLIAMO LA MATRICE A (versione scalare, poi va "duplicata")
    % ------------------------------------------------------
    % gradi di liberta' globali:
    %
    % vertice i  ----> i
    %
    % elemento l -----> nver+l
    %
    %
    % array dei dof globali
    %
    dofg = [v1 v2 v3 nver+iele];   % 1a componente
    dofgg = [dofg , dofg + (nver+nele)];  % aggiungo anche 2a componente    
    %   
    A(dofgg,dofgg) = A(dofgg,dofgg) + KE;     % assemblo
    %
    % termine noto
    %
    fE1 = zeros(4,1);   % carico per la prima componente
    fE2 = zeros(4,1);   % carico per la seconda componente
    %
    for i=1:4
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
   dofp = [v1 v2 v3];    % GdL globali della pressione (indici dei vertici)
   B(dofgg,dofp) = B(dofgg,dofp) + BE;     
end

b = [b1;b2];    % combino le due componenti del carico in singolo vettore


% ---------------------------------
% sistemiamo le condizioni al bordo
% ---------------------------------

internV = setdiff([1:1:nver],boundary); % indici dei vertici INTERNI

% Ora considero che i GDL sono per prima e poi seconda componente:
NL  = [internV,(nver+nele) + internV];           % tutti i GdL interni

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
Kh = [Ah,Bh;Bh',zeros(nver)];
fh = [fh;zeros(nver,1)];


% -----------------------------------------------------
% SE le condizioni al bordo sulla velocita' sono tutte
% ... Dirichlet omogenee, si deve imporre condizione 
% ... media nulla sulla pressione. Lo facciamo "con moltiplicatore".
% ALTRIMENTI commentare le righe sottostanti.
% -----------------------------------------------------

N = length(NL);
Kh = [Kh,[zeros(N,1);vettarea];zeros(1,N),vettarea',0];
fh = [fh;0];


% -------------------------------
% Risolviamo il sistema lineare
% -------------------------------
uh = zeros(2*(nver+nele),1);
solh = Kh\fh;
uh(NL) = solh(1:N);            % estraggo le velocita'
uh1 = uh(1:length(uh)/2);      % estraggo la prima componente
uh2 = uh(length(uh)/2+1:end);  % estraggo la seconda componente
ph  = solh(N+1:N+nver);        % estraggo le pressioni


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
phihq = zeros(4,Nq);
%
for i=1:4
    for q=1:Nq
        phihq(i,q) = phih2(i,xhq(q),yhq(q));
        [gx,gy] = gradhphih2(i,xhq(q),yhq(q));
        gphihqx(i,q) = gx;
        gphihqy(i,q) = gy;
    end
end
%
errL2sq = 0;   % errore L2 (velocit?)
errH1sq = 0;   % errore H1 (velocit?)
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
    % elemento l -----> nver+l
    %
    % array dei dof globali del triangolo T
    %
    dofg = [v1 v2 v3 nver+nele];
    %
    uT1 = uh1(dofg);   % prima componente velocit? (locale)
    uT2 = uh2(dofg);   % seconda componente velocit? (locale)
    pT = ph([v1 v2 v3]); % pressione (locale)
    %
    sq  = 0;
    sqH = 0;
    sqP = 0;
    %
    for q=1:Nq
        %
        % calcolo la sommatoria sulle phi_i
        %
        tmp1  = 0;       % valori della funzione uh (prima componente)
        tmp2  = 0;       % (seconda componente)
        tmpH1 = [0;0];   % valori dei gradienti di uh (prima componente)
        tmpH2 = [0;0];   % (seconda componente)
        tmpP = 0;        % valori della pressione ph
        %
        for i=1:4
            tmp1 = tmp1 + uT1(i)*phihq(i,q);
            tmp2 = tmp2 + uT2(i)*phihq(i,q);
            tmpH1= tmpH1 + uT1(i)*(JFIT*[gphihqx(i,q);gphihqy(i,q)]); 
            tmpH2= tmpH2 + uT2(i)*(JFIT*[gphihqx(i,q);gphihqy(i,q)]);
        end
        %
        tmpP = pT(1)*phihq(1,q) + pT(2)*phihq(2,q) + pT(3)*phihq(3,q);
        %
        pos = JF*[xhq(q);yhq(q)]+[x1;y1];  
        %
        xq = pos(1);         % x-coordinata del gauss-nodo q-esimo
        yq = pos(2);         % y-coordinata del gauss-nodo q-esimo
        %
        sq = sq + (uexact(xq,yq,1)-tmp1)^2*whq(q) ...
                + (uexact(xq,yq,2)-tmp2)^2*whq(q);     
        sqH = sqH + norm(uexactG(xq,yq,1) - tmpH1)^2*whq(q) ...
                  + norm(uexactG(xq,yq,2) - tmpH2)^2*whq(q);
        sqP = sqP + (pexact(xq,yq) - tmpP)^2*whq(q);       
    end
    %
    sq = sq*2*area;
    sqH = sqH*2*area;
    sqP = sqP*2*area;
    %
    errL2sq = errL2sq + sq;
    errH1sq = errH1sq + sqH;
    preL2sq = preL2sq + sqP;
    %
end
%
errL2 = sqrt(errL2sq);  % (comprende ambo le componenti)
errH1 = sqrt(errH1sq);  % (comprende ambo le componenti) 
errPR = sqrt(preL2sq);  % (errore pressioni)
%


