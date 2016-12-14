%
function [uh1,uh2] = stokesP2(xv,yv,vertices,edges,endpoints,boundary,boundedges);    %boundedges= nodi dei punti medi
% ---------------------------
% risolve il problema di Stokes con polinomi di grado 2
% plotta i grafici della soluzione (vettoriale) e con la funzione quiver

% Ap2 matrice di stiffness
% ---------------------------

nver = length(xv);          % calcolo numero vertici
nedge = size(endpoints,1);  % calcolo numero lati
nele = size(vertices,1);    % calcolo numero elementi

internV = setdiff([1:1:nver],boundary); % indici dei vertici INTERNI, in caso anche alcuni di bordo senza vincoli
internE = setdiff([1:1:nedge],boundedges); % indici dei lati INTERNI
NL = [internV,nver+internE];     % gradi di libertà interni

% -------------------------------
% CALCOLO DELLA MATRICE A
% -------------------------------

%
% FORMULA DI QUADRATURA
%
fdq = 'degree=4';
%
[xhq,yhq,whq]=quadratura(fdq);
%
Nq = length(xhq);
%
%
% (xhq,yhq) = nodi di quadratura su T cappello
% whq = relativi pesi
%
% funzioni di base calcolate nei nodi di
% quadratura dell'elemento di riferimento (che restano uguali in tutto il processo di risoluzione)
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
% nei nodi di quadratura dell'elemento di riferimento (anche i gradienti rimangono uguali)
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
% ASSEMBLIAMO LA MATRICE DEI COEFFICIENTI (STIFFNESS) -blocco "base" del sistema stazionario con 
% polinomi di grado due da posizionare a nord ovest e sud est
%
Ap2 = sparse(nver+nedge,nver+nedge);

%
for iele=1:nele     %ciclo su tutti gli elementi della mesh
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
    KE = zeros(6,6);    %matrice grad grad locale
    
    %
    for i=1:6       %sfrutto la simmetria delle due matrici
        for j=1:i-1
            KE(i,j) = KE(j,i);      %matrice stifness locale
            
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
                xq = tmp(1);    % x-coordinata del nodo di quad.
                yq = tmp(2);    % y-coordinata del nodo di quad.
                %
                tmpK = dot(...
                    (JFIT*[gphihqx(j,q);gphihqy(j,q)]),...
                    (JFIT*[gphihqx(i,q);gphihqy(i,q)])...
                    );
                %
                %
                KE(i,j) = KE(i,j) + c(xq,yq)*tmpK*whq(q);   %sommo tutti i contributi delle matrici locali
            end
            %
            KE(i,j) = 2*area*KE(i,j);
            %
        end
    end
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
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];  %recupera gli indici dei gradi di libertà guardando la situazione della mesh globale
    %
    Ap2(dofg,dofg) = Ap2(dofg,dofg) + KE;   %sto sommando i contributi della matrice locale alla matrice grande
    
end

Ap2 = Ap2(NL,NL);   % ritaglio la matrice nei dofs interni


%
% costruisco la matrice "grande" A con i blocchi costruiti Ap2
%

A = zeros( 2*length(NL), 2*length(NL));
A(1:length(NL),1:length(NL)) = Ap2;
A((length(NL)+1):2*length(NL),(length(NL)+1):2*length(NL)) = Ap2;



%
% calcolo il carico del problema componente per componente
%
%aggiungo il parametro k alla function P2load e a carico in modo da calcolare i 
%vettori f1 ed f2 con un'unica function e non 
%con la stessa differente per una sola istruzione!

b1 = P2load(xv,yv,vertices,edges,boundary,boundedges,endpoints,1);    
b2 = P2load(xv,yv,vertices,edges,boundary,boundedges,endpoints,2);

b=[b1; b2];


uh = A\b;

%----------------
% La uh e' un vettore con i valori solo sui dofs interni;
% ricostruisco la vera uh aggiungendo valore zero ai dofs di bordo
% ----------------


%separo le componenti di uh 
uh1 = uh(1:length(NL));
uh2 = uh(length(NL)+1:2*length(NL));


% ----------------
% Le uh1 e uh2 sono vettori con i valori solo sui dofs interni;
% ricostruisco ciascuna componente di uh aggiungendo valore zero ai dofs di bordo
% ----------------

uhh1 = uh1;
uh1 = zeros(nver+nedge,1);
uh1 (NL) = uhh1;

uhh2 = uh2;
uh2 = zeros(nver+nedge,1);
uh2 (NL) = uhh2;

uh = [uh1 uh2]; % costruisco una matrice unica con prima colonna contenente la prima componente uh1
% e nella seconda colonna la seconda componente in modo da usare un unico
% testo nel programma



% ---------------------------------
% CALCOLO DEGLI ERRORI
% ---------------------------------

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

% componente 1
errL2sq1 = 0;   % errore L2 
errH1sq1 = 0;   % errore H1
uH1sq1 = 0;
%
% componente 2 
errL2sq2 = 0;   % errore L2 
errH1sq2 = 0;   % errore H1
uH1sq2 = 0;



for j=1:2   % ciclo sulla colonna j-esima, cioè sulla componente j-esima della soluzione uh
    %calcolando i due errori separatamente
    
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
        uT = uh(dofg,j);
        %
        sq  = 0;
        sqH = 0;
        uH = 0;   % serve per norma di sol esatta (ai fini di fare errore relativo)
        %
        for q=1:Nq
            %
            % calcolo la sommatoria sulle phi_i
            %
            tmp  = 0;       % valori della funzione uh
            tmpH = [0;0];   % valori dei gradienti di uh
            %
            for i=1:6
                tmp = tmp + uT(i)*phihq(i,q);   % è la uh nei nodi di Gauss
                tmpH= tmpH + uT(i)*(JFIT*[gphihqx(i,q);gphihqy(i,q)]);
            end
            %
            tmp2 = JF*[xhq(q);yhq(q)]+[x1;y1];
            %
            xq = tmp2(1);
            yq = tmp2(2);
            %
            sq = sq + (uexact(xq,yq,j)-tmp)^2*whq(q);
            sqH = sqH + norm(uexactG(xq,yq,j) - tmpH)^2*whq(q);
            uH = uH + norm(uexactG(xq,yq,j))^2*whq(q);
            %
        end
        %
        sq = sq*2*area;
        sqH = sqH*2*area;
        uH = uH*2*area;
        %
        
        if j==1     % considero la prima componente e memorizzo l'errore calcolato
            
            errL2sq1 = errL2sq1 + sq;
            errH1sq1 = errH1sq1 + sqH;
            uH1sq1 = uH1sq1 + uH;
            
        else        %considero la seconda componente e memorizzo l'errore calcolato
            
            errL2sq2 = errL2sq2 + sq;
            errH1sq2 = errH1sq2 + sqH;
            uH1sq2 = uH1sq2 + uH;
            
        end
        %
    end
    %
    
end

disp('errore (assoluto) in norma L2 (rispettivamente prima e seconda componente)')
errL2 = sqrt(errL2sq1)
errL2 = sqrt(errL2sq2)
%
disp('errore (relativo) in norma H1 (rispettivamente prima e seconda)')
errH1 = sqrt(errH1sq1/uH1sq1)
errH1 = sqrt(errH1sq2/uH1sq2)
%

% -----------
% PLOT DELLA SOLUZIONE uh (basato sui valori ai vertici)
% -----------

% plot fatto con ciclo for per evitare doppio testo

for i=1:2
    
    figure(i)
    
    uhi = uh(:,i);      % sto memorizzando la colonna i-esima della matrice uh
                        % (cioè sto considerando la componente i-esima della soluzione uh)
    
    for k=1:size(vertices,1)
        hold on;
        index=(vertices(k,1:3))';
        % riga sotto: opzione per funzione nota ai vertici
        uu_tmp=[uhi(index);uhi(index(1))];
        % riga sotto: opzione per funzione costante a tratti
        % p_tmp = ph(k)*ones(length(index)+1,1);
        vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
        fill3(vert_temp(:,1),vert_temp(:,2),uu_tmp,uu_tmp);
    end
    view(3)
    grid on
    colorbar
    hold off
    
end


%-------
%plotto la soluzione uh in forma vettoriale tramite comando quiver
%-------

v = zeros(nver,1);      %sto inizializzando un vettore con tutti gli indici di numerazione dei vari vertici

for iele = 1:nele
    %
    %ciclo su tutti gli elementi della mesh ed estraggo gli indici dei
    %verici creando un array in base a cui tagliare le componenti della soluzione
    %uh
    %
    v1 = vertices(iele,1);
    v2 = vertices(iele,2);
    v3 = vertices(iele,3);
    
    x1 = xv(v1);
    y1 = yv(v1);
    %
    x2 = xv(v2);
    y2 = yv(v2);
    %
    x3 = xv(v3);
    y3 = yv(v3);
    %
    
    v (v1) = v1;
    v (v2) = v2;
    v (v3) = v3;
    
    
end


figure(3)

quiver( xv, yv, uh1(v), uh2(v));    
        
        
       
