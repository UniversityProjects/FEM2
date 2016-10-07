function [A,M] = P1stiffmass(xv,yv,vertices,edges,boundary,endpoints);
%
% Modifica del file dei P1 che hanno gia da Ale Russo, che lavora
% ... come spiegato nel file readme di questa cartella. Calcola le
% ... matrici stiffness e massa dello stazionario (diffusione).
%
area = 0;
%
nver = length(xv); 
nele = size(vertices,1);
%
A = sparse(nver,nver);
M = sparse(nver,nver);
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
    % matrici elementari
    %
    Aloc = 1/(4*T)*[dot(e1,e1) dot(e2,e1) dot(e3,e1) 
                  dot(e1,e2) dot(e2,e2) dot(e3,e2)
                  dot(e1,e3) dot(e2,e3) dot(e3,e3)];
    %    
    % A = c(xb,yb)*A;  % c scalare davanti al grad-grad (OPZIONALE)
    %
    Mloc = (T/12)*ones(3) + (T/12)*eye(3);  % massa locale (usare pti medi)
    % Mloc = T/9*ones(3,3);  
    %    
    % assemblo le matrici globali
    %
    A([v1 v2 v3],[v1 v2 v3]) = ...
        A([v1 v2 v3],[v1 v2 v3]) + Aloc;
    M([v1 v2 v3],[v1 v2 v3]) = ...
        M([v1 v2 v3],[v1 v2 v3]) + Mloc;
end
%
% Adesso sistemiamo le condizioni al bordo (omogenee Dirichlet)
% ... e ci costruiamo le vere A e M 
%
intern = setdiff([1:1:nver],boundary); % indici dei nodi INTERNI
A = A(intern,intern);
M = M(intern,intern);







    
    
    
    











        
        
        




















    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





























