function bh = P1load(xv,yv,vertices,edges,boundary,endpoints,t);
%
% Modifica del file dei P1 che hanno gia da Ale Russo, che lavora
% ... come spiegato nel file readme di questa cartella. Calcola il
% ... vettore termine noto (dato tramite carico.m) al tempo t.
%
area = 0;
%
nver = length(xv); % numero di vertici
nele = size(vertices,1);
%
bh = zeros(nver,1);      
%
% queste sono matrici e vettori SENZA bc
%
% facciamo un loop su tutti gli elementi
%
for iele=1:nele
    %
    % recuperiamo le informazioni dell'elemento iele-esimo
    %
    % vertici
    v1 = vertices(iele,1);
    v2 = vertices(iele,2);
    v3 = vertices(iele,3);
    %
    % coordinate dei vertici
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
    % vettori dei lati
    e1 = [x3-x2,y3-y2];
    e2 = [x1-x3,y1-y3];
    e3 = [x2-x1,y2-y1];
    %
    % area
    T = 0.5* det([1  1  1
                 x1 x2 x3
                 y1 y2 y3]);
    %
    % termine noto elementare
    %
    % trapezi 
    fhT_t = T/3 * [carico(x1,y1,t) carico(x2,y2,t) carico(x3,y3,t)]';
    %
    % baricentro
    fhT_b = T/3 * carico(xb,yb,t) * [1 1 1]';
    %
    % fdq esatta per polinomi di grado 2
    %
    % cominciamo col definire i punti medi
    % dei lati
    %
    M1x = (x2+x3)/2;
    M1y = (y2+y3)/2;
    %
    M2x = (x3+x1)/2;
    M2y = (y3+y1)/2;
    %
    M3x = (x1+x2)/2;
    M3y = (y1+y2)/2;
    %
    fhT_2 = T/6*[carico(M2x,M2y,t)+carico(M3x,M3y,t)
                 carico(M1x,M1y,t)+carico(M3x,M3y,t)
                 carico(M1x,M1y,t)+carico(M2x,M2y,t)];
    %
    % assembliamo il termine noto globale
    %
    bh([v1 v2 v3]) = bh([v1 v2 v3]) + fhT_2; 
    %
end

% adesso ritagliamo il carico ai nodi interni 
% (cond al bordo: Dirichlet omogenee)

intern = setdiff([1:1:nver],boundary); % indici dei nodi INTERNI
bh = bh(intern);



    
    
    
    











        
        
        




















    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





























