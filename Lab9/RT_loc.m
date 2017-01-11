function [K,P,f,M,E]=RT_loc(Xv,Yv);
%---------------
% data matrice X che determina (come sotto) la geometria del triangolo,
% ... questo algo calcola ma matrice 3x3 di stiffness K per gli elementi RT
% --------------
% INPUT:
% Vettore Xv riga con le x-coordinate dei tre vertici. E analogo vettore Yv
% ... per le y-coord. I vertici sono assunti ordinati in senso anti-orario
%----------------
% OUTPUT:
% matrice 3x3 K di massa dei flussi
% matrice 3x1 P di pressioni contro divergenza dei flussi
% matrice M di trasformazione da "scrittura dofs" a "scrittura [a,b,c]"
% scalare E area del elemento
% f scalare che rappresenta il carico locale
%----------------
% NOTE:
% Non usa trasf. di Piola, conti fatti direttamente sul elemento fisico.
% I gradi di liberta' sono i soliti valori normali su ogni lato.
% Si puo scrivere ognuna di tali funzioni dello spazio RT anche nella
% ... forma [a,b]^T + c [x,y]^T ...
% ... con x,y coordinate centrate nel baricentro del elemento
%----------------
X = [Xv;Yv];
% (X(:,2)-X(:,1))/norm(X(:,2)-X(:,1)) tangente al primo lato 
% costruisco N matrice delle normali 
% (i-esima colonna) tiene la normale al i-esimo lato
L = [norm(X(:,3)-X(:,2));norm(X(:,1)-X(:,3));norm(X(:,2)-X(:,1))];
% L e' un vettore con la lunghezza dei tre lati
N=[(X(:,3)-X(:,2))/L(1), ...
    (X(:,1)-X(:,3))/L(2),...
    (X(:,2)-X(:,1))/L(3)];
N = [0,1;-1,0]*N;
% ---
% calcolo M, matrice che trasforma il vettore dei gradi di liberta' 
% ... nei valori a,b,c associati alla scrittura di cui sopra
P = [X(:,2)+X(:,3),X(:,3)+X(:,1),X(:,1)+X(:,2)]/2; % matrice dei punti medi
b = ((sum(X'))')/3;     % baricentro del elemento
TM = P - [b,b,b];   % scalo via il baricentro dai punti medi 
M = [N' , diag((TM')*N)];
M = inv(M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = cross([(X(:,2)-X(:,1));0],[(X(:,3)-X(:,1));0]);  
E = E(3)/2;             % area del elemento
II = (E/3)*sum(sum((TM).^2)); 
% II e' l'integrale di x^2+y^2 sul elemento (scalati col baricentro) 
% ... calcolato usando la formula dei punti medi del triangolo
A = [E*eye(2),zeros(2,1);zeros(1,2),II];
% A e' la matrice di stiffness calcolata rispetto ai coeffs [a,b,c]
K = M'*A*M;

% Calcolo vettore 3x1 P (conto immediato andando per parti)
P = L;   

% Calcolo carico f (il segno meno segue dalla formulazione variazionale)
f = - E*carico(b);    % integrale tramite valutazione nel baricentro






