%
function [errL2, errH1] = P2_Vect_vol (...
                                       xv, yv, ...
                                       vertices, ...
                                       edges, ...
                                       endpoints, ...
                                       boundary, ...
                                       boundedges, ...
                                       lambda, ...
                                       out, plot ...
                                       )
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
B = sparse(2*(nver+nedge),2*(nver+nedge));
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
    % baricentro del triangolo
    xb = (x1+x2+x3)/3;
    yb = (y1+y2+y3)/3; 
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
    area = 0.5*det(JF);
    %
    KE = zeros(6,6);    % matrici locali per la A (vedi NOTA in cima)   
    BE = zeros(12,12);   % matrici locali per la B
    %
    for i=1:6
%         for j=1:i-1
%             KE(i,j) = KE(j,i);
%         end  
        for j=1:6
            for q=1:Nq
                %
                % calcoliamo l'immagine (xq,yq) sul triangolo
                % corrente del nodo di quadratura
                % (xhq(q),yhq(q)) che sta 
                % sull'elemento di rif.
                tmp = dot(... % Grad:Grad Term
                      (JFIT*[gphihqx(j,q);gphihqy(j,q)]),...
                      (JFIT*[gphihqx(i,q);gphihqy(i,q)])...
                      );
               %
                KE(i,j) = KE(i,j) + tmp*whq(q);
               %  KE(i,j) = KE(i,j) + c(xq,yq)*tmp*whq(q);            
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
        for i=1:12
            if j < 6.5         % prime sei funzioni di base
                [gx, gy] = gradhphih2 (j, xb, yb);
                tmpJ = JFIT(1,1:2)*[gx; gy];
            elseif j > 6.5   % seconde sei funzioni di base
                [gx, gy] = gradhphih2 (j-6, xb, yb);
                tmpJ = JFIT(2,1:2)*[gx; gy];
            end
            if i < 6.5         % prime sei funzioni di base
                [gx, gy] = gradhphih2 (i, xb, yb);
                tmpI = JFIT(1,1:2)*[gx; gy];
            elseif i > 6.5   % seconde sei funzioni di base
                [gx, gy] = gradhphih2 (i-6, xb, yb);
                tmpI = JFIT(2,1:2)*[gx; gy];
            end
            vol = lambda*tmpJ*tmpI;
            BE(i, j) = area*vol;   
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
    dofg = [v1 v2 v3 nver+l1 nver+l2 nver+l3];   % prima componente
    dofgg = [dofg , dofg + (nver+nedge)];     % aggiungo seconda componente
    %
    A(dofgg,dofgg) = A(dofgg,dofgg) + KE;
    %
    % ------------------------------------------------------
   % ASSEMBLIAMO LA MATRICE B
   % ------------------------------------------------------
   B(dofgg,dofgg) = B(dofgg,dofgg) + BE;
end

% ---------------------------------
% CALCOLO IL TERMINE DI CARICO 
% Le due componenti fatte separatamente usando P2load.m.
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


Ah = A(NL,NL);
Bh = B(NL,NL);
fh = b(NL);

% disp(['sizeA = ', num2str(size(A))]);
% disp(['sizeAh = ', num2str(size(Ah))]);
% disp(['sizB = ', num2str(size(B))]);
% disp(['sizeBh = ', num2str(size(Bh))]);
% disp(['sizeb = ', num2str(size(b))]);
% disp(['sizefh = ', num2str(size(fh))]);



%
clear A
clear B
clear b
% -----------------------------------------------------
% Costruisco la matrice totale (parti A e B) e il carico
%
Kh = Ah + Bh;
% disp(['sizeKh = ', num2str(size(Kh))]);
% disp(['sizefh = ', num2str(size(fh))]);

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
% Adesso calcoliamo l'errore in L2 e H1 rispetto a ue
%
[errL2, errH1] = err(uh, 'degree=5', xv, yv, vertices, edges, endpoints, out);

%

% ----------------------------------
% PLOT DELLA SOLUZIONE uh COMPONENTE PER COMPONENTE
% (basato sui valori ai vertici)
% ----------------------------------
if (strcmp(plot,'yes'))
figure()
for k=1:size(vertices,1)
    hold on;
    index=(vertices(k,1:3))';
    % riga sotto: opzione per funzione nota ai vertici
    uu_tmp=[uh1(index);uh1(index(1))];  % prima componente
    % uu_tmp=[uh2(index);uh2(index(1))];  % seconda componente
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


% ----------------------------------
% PLOT DELLA SOLUZIONE u ESATTA COMPONENTE PER COMPONENTE
% (basato sui valori ai vertici)
% ----------------------------------
if (strcmp(plot,'yes'))
figure()
fmesh(@(x,y) 0.25*((x^2 - 1)^2)*(y^2 - 1)*y, [-1 1 -1 1]);
title('Exact Solution');
view(3)
grid on
colorbar
hold off    
end

        
% ----------------------------------
% PLOT DELLA SOLUZIONE uh COME FRECCE AI VERTICI
% ----------------------------------        
% 
if (strcmp(plot,'yes'))
    figure()
    quiver(xv,yv,uh1(1:nver),uh2(1:nver));
end
%
%




