function [uh] = stationary(xv,yv,vertices,edges,boundary,endpoints);
%
% Algo PRELIMINARE che risolve il pbl della diffusione stazionario
% Elementi P1 continui (vedi readme file per link con corso A.Russo)
% ---------------------------
% Il carico (dipendente dal tempo) e' dato nel file carico.m
% A matrice di stiffness
% M matrice di massa
% ---------------------------
nver = length(xv); 
intern = setdiff([1:1:nver],boundary); % indici dei vertici INTERNI
nint=length(intern);                   % numero vertici interni
%
%--------------------------------
% costruisco matrici di stiffness e massa
[A,M] = P1stiffmass(xv,yv,vertices,edges,boundary,endpoints);
% ------------------------------

alfa=1;  % fattore davanti alla diffusione
gamma=0; % fattore davanti alla reazione

MAT = alfa*A + gamma*M;  % calcolo matrice del sistema

rhs = P1loadS(xv,yv,vertices,edges,boundary,endpoints); % calcolo carico
%----------
uh = MAT\rhs;     % risoluzione sistema lineare
%----------

% --------
% La uh e' un vettore con i valori solo sui nodi interni;
% ricostruisco la vera uh aggiungendo valore zero ai nodi di bordo
uhh=uh;     
uh=zeros(nver,1);
uh(intern)=uhh;
% --------

% --------
% ERRORE IN NORMA DEL MASSIMO AI NODI
% --------
uI = exsol(xv,yv);
disp('Errore massimo ai nodi')
err = max(abs(uI-uh))/max(abs(uI))

% --------
% ERRORE L2 DISCRETA (valutata ai vertici,valido per mesh quasi-uniformi)
% --------
% disp('Errore in norma L2 discreta ai nodi')
% errL2 = norm(uI-uh)/norm(uI)

% --------
% ERRORE NORME DISCRETE con A,M 
% --------
uI=uI(intern);
disp('Errore in norme L2 ed H1 discrete')
errM = sqrt( ((uI-uhh)'*M*(uI-uhh))/(uI'*M*uI) )
errA = sqrt( ((uI-uhh)'*A*(uI-uhh))/(uI'*A*uI) )


% --------
% ERRORE SEMI-NORMA H1 
% --------

errH = 0;
unorm = 0;
nele = size(vertices,1);
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
    % valori della soluzione discreta ai vertici
    uh1 = uh(v1);
    uh2 = uh(v2);
    uh3 = uh(v3);
    % recupero il valore del gradiente di uh sul elemento
    QMAT = [e1;e2];
    Qvec = [uh3-uh2;uh1-uh3];
    gradh = QMAT\Qvec; 
    % calcolo errore sul elemento e lo aggiungo al totale
    errH = errH + T*norm(exsolG(xb,yb) - gradh)^2;  
    unorm = unorm + T*norm(exsolG(xb,yb))^2; 
end
disp('errore in norma H1')
errH = sqrt(errH/unorm)
    

% -----------
% PLOT DELLA SOLUZIONE uh
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







    
    