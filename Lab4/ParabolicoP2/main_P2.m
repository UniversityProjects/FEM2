%
function [uh] = main_P2(xv,yv,vertices,edges,endpoints,boundary,boundedges);
%
% Risolve il problema della diffusione dipendente dal tempo
% Elementi P2 continui
% ---------------------------
% Il carico (dipendente dal tempo) e' dato nel file carico.m
% Il dato iniziale e' dato nel file initial.m
% A matrice di stiffness
% M matrice di massa
% ---------------------------
T = 1;   % tempo finale, il tempo iniziale e' sempre zero
N =5; % numero di intervalli temporali
dt= T/N; % delta t
th = 1; % valore theta del theta-metodo 

nver = length(xv);          % calcolo numero vertici
nedge = size(endpoints,1);  % calcolo numero lati
nele = size(vertices,1);    % calcolo numero elementi

internV = setdiff([1:1:nver],boundary); % indici dei vertici INTERNI
internE = setdiff([1:1:nedge],boundedges); % indici dei lati INTERNI
NL = [internV,nver+internE];     % gradi di libertà interni

% ------------------------
% costruisco dato iniziale
% ------------------------
uh=zeros(length(NL),1);
%
% scorre i gradi di libertà interni e calcola dato iniziale (interpolato) 
%
for k=1:length(internV)            % CICLO SUI VERTICI INTERNI
    uh(k) = initial(xv(internV(k)),yv(internV(k)));   
end
%
for k=1:length(internE)            % CICLO SUI LATI INTERNI
    v1 = endpoints(internE(k),1);  % indice primo endpoint  
    v2 = endpoints(internE(k),2);  % indice secondo endpoint
    xm = (xv(v1)+xv(v2))/2;        % x-coord del pto medio del lato
    ym = (yv(v1)+yv(v2))/2;        % y-coord del pto medio del lato
    uh(length(internV)+k) = initial(xm,ym);   
end
%


% -------------------------------
% PRE-CALCOLO DELLE MATRICI A e M
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
% ASSEMBLIAMO LA MATRICE DEI COEFFICIENTI (STIFFNESS E MASSA)
% 
A = sparse(nver+nedge,nver+nedge);
M = sparse(nver+nedge,nver+nedge);

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
    ME = zeros(6,6);
    %
    for i=1:6
        for j=1:i-1
            KE(i,j) = KE(j,i);
            ME(i,j) = ME(j,i);
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
                tmpM = phihq(j,q)*phihq(i,q);
                %
                KE(i,j) = KE(i,j) + c(xq,yq)*tmpK*whq(q);
                ME(i,j) = ME(i,j) + tmpM*whq(q);
            end
            %
            KE(i,j) = 2*area*KE(i,j);
            ME(i,j) = 2*area*ME(i,j);
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
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];
    %
    A(dofg,dofg) = A(dofg,dofg) + KE;
    M(dofg,dofg) = M(dofg,dofg) + ME;
end

A = A(NL,NL);   % ritaglio la matrice nei dofs interni
M = M(NL,NL);   % ritaglio la matrice nei dofs interni

% ---------------------------------
% ADESSO INIZIAMO LE ITERAZIONI 
% ---------------------------------

for j=0:N-1
    t = (j + th)*dt;
    % calcolo il vettore di carico al tempo t = (j + theta)*dt
    b = P2load(xv,yv,vertices,edges,boundary,boundedges,endpoints,t);
    % ---
    rhs = ( M + (th-1)*dt*A )*uh + dt*b;  % membro destro del theta-metodo
    MAT = (M + th*dt*A);    % matrice del theta metodo
    %----------    
    uh = MAT\rhs;     % risoluzione sistema lineare al passo
    %----------
end


% ----------------
% La uh e' un vettore con i valori solo sui dofs interni;
% ricostruisco la vera uh aggiungendo valore zero ai dofs di bordo
% ----------------
uhh=uh;     
uh = zeros(nver+nedge,1);
uh(NL)=uhh;


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
errL2sq = 0;   % errore L2
errH1sq = 0;   % errore H1
uH1sq = 0;
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
    uT = uh(dofg);
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
            tmp = tmp + uT(i)*phihq(i,q);
            tmpH= tmpH + uT(i)*(JFIT*[gphihqx(i,q);gphihqy(i,q)]); 
        end
        %
        tmp2 = JF*[xhq(q);yhq(q)]+[x1;y1];
        %
        xq = tmp2(1);
        yq = tmp2(2);
        %
        sq = sq + (exsol(xq,yq)-tmp)^2*whq(q);
        sqH = sqH + norm(exsolG(xq,yq) - tmpH)^2*whq(q);
        uH = uH + norm(exsolG(xq,yq))^2*whq(q);
        %
    end
    %
    sq = sq*2*area;
    sqH = sqH*2*area;
    uH = uH*2*area;
    %
    errL2sq = errL2sq + sq;
    errH1sq = errH1sq + sqH;
    uH1sq = uH1sq + uH;
    %
end
%
% disp('errore (assoluto) in norma L2')
% errL2 = sqrt(errL2sq)
disp('errore (relativo) in norma H1')
errH1 = sqrt(errH1sq/uH1sq)
%

% -----------
% PLOT DELLA SOLUZIONE uh (basato sui valori ai vertici)
% -----------
figure(1)
for k=1:size(vertices,1)
    hold on;
    index=(vertices(k,1:3))';
    % riga sotto: opzione per funzione nota ai vertici
    uu_tmp=[uh(index);uh(index(1))];
    % riga sotto: opzione per funzione costante a tratti
    % p_tmp = ph(k)*ones(length(index)+1,1); 
    vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
    fill3(vert_temp(:,1),vert_temp(:,2),uu_tmp,uu_tmp); 
end
view(3)
grid on
colorbar
hold off        
        
        
       
