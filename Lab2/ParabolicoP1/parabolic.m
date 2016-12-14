function [uh] = parabolic(xv,yv,vertices,edges,boundary,endpoints);
%
% Algo PRELIMINARE che risolve il pbl della diffusione stazionario
% Elementi P1 continui (vedi readme file per link con corso A.Russo)
% ---------------------------
% Il carico (dipendente dal tempo) e' dato nel file carico.m
% A matrice di stiffness
% M matrice di massa
% ---------------------------

clc

%--------------------------------
% Mesh Generation
%--------------------------------
hmax = 0.1;
[xv,yv,vertices,edges,endpoints,boundary,boundedges] = meshgen(hmax);
disp(' --- Mesh Generation --- ');
%
%--------------------------------
nver = length(xv); 
intern = setdiff([1:1:nver],boundary); % indici dei vertici INTERNI
nint=length(intern);                   % numero vertici interni
%
% ------------------------------
% Time Discretization
% ------------------------------
disp(' --- Time Discretization --- ');
t_size = 20;
t_0 = 0;
t_max = 1;
tau = t_size / t_max;
[t] = time(t_size, t_0, t_max);
%
%--------------------------------
% Costruisco Le Matrici
disp(' --- Matrrices Building --- ');
% Costruisco le matrici A e M
[A,M] = P1stiffmass(xv,yv,vertices,edges,boundary,endpoints);
% Parametro theta per il metodo di risoluzione
% theta = 0; % Eulero In Avanti (Esplicito)
% theta = 0.5; % Punto Medio
theta = 1; % Eulero All'Indietro (Implicito)
%
% Costruisco la matrice M + tau*theta*A
l_matrix = M + tau*theta*A;
%disp('size_l = '),disp(size(l_matrix));
% Costruisco la matrice M - tau*(1-theta)*A
r_matrix = M - tau*(1-theta)*A;
%disp('size_r = '),disp(size(r_matrix));
%
% ------------------------------
% Calcolo Soluzione Sistema
% ------------------------------
disp(' --- System Solution --- ');
% Definition of the solution array
uht = zeros(t_size,nint);
%disp('size_uht = '),disp(size(uht));
% Add initial condition on the internal nodes
for j=1:nint
    uht(1,j) = initial(xv(intern(j)),yv(intern(j)));
end
%disp('size_uht(i,:) = '),disp(size(uht(1,:)));
for i=2:t_size
    time_load = (1-theta)*t(i-1) + theta*t(i);
    load = P1loadS(xv,yv,vertices,edges,boundary,endpoints,time_load); % calcolo carico
    %disp('size_load = '),disp(size(load));
    %disp('size_prod = '),disp(size(r_matrix*uht(i-1,:)'));
    uht(i,:) = l_matrix \ ( r_matrix*uht(i-1,:)' + tau*load);
end
%
% ------------------------------
% Aggiunta Condizioni Al Bordo
% ------------------------------
% uht è solo aui nodi interni, aggiungo le condizioni sui nodi di bordo
disp(' --- Border Condition --- ');
uhtt = uht;   
uht = zeros(t_size,nver);
uht(:, intern) = uhtt;
%-----------------------------
%
%-----------------------------
% ERRORE IN NORMA DEL MASSIMO AI NODI
%-----------------------------
disp(' --- Errors --- ');
uI = exsol(xv,yv);
err = max(abs(uI-uht(t_size,:)'))/max(abs(uI));
disp('Errore massimo ai nodi');
disp(['  err = ', num2str(err)]);
%
%-----------------------------
% ERRORE L2 DISCRETA (valutata ai vertici,valido per mesh quasi-uniformi)
% --------
errL2 = norm(uI-uht(t_size,:)')/norm(uI);
disp('Errore in norma L2 discreta ai nodi');
disp(['  errL2 = ', num2str(errL2)]);
%
%-----------------------------
% ERRORE NORME DISCRETE con A,M 
%-----------------------------
uI=uI(intern);
disp('Errore in norme L2 ed H1 discrete')
errM = sqrt( ((uI-uhtt(t_size,:)')'*M*(uI-uhtt(t_size,:)'))/(uI'*M*uI) );
disp(['  ErrM = ', num2str(errM)]);
errA = sqrt( ((uI-uhtt(t_size,:)')'*A*(uI-uhtt(t_size,:)'))/(uI'*A*uI) );
disp(['  ErrA = ', num2str(errA)]);
%
%-----------------------------
% ERRORE SEMI-NORMA H1 
%-----------------------------
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
    uht1 = uht(t_size,v1);
    uht2 = uht(t_size,v2);
    uht3 = uht(t_size,v3);
    % recupero il valore del gradiente di uh sul elemento
    QMAT = [e1;e2];
    Qvec = [uht3-uht2;uht1-uht3];
    gradh = QMAT\Qvec; 
    % calcolo errore sul elemento e lo aggiungo al totale
    errH = errH + T*norm(exsolG(xb,yb) - gradh)^2;  
    unorm = unorm + T*norm(exsolG(xb,yb))^2; 
end
errH = sqrt(errH/unorm);
disp('Errore in norma H1');
disp(['  errH1 = ', num2str(errH)]);
%
%-----------------------------
% PLOT DELLA SOLUZIONE uht
%-----------------------------
disp(' --- Plot uht ---');
figure(2)
for k=1:size(vertices,1)
    hold on;
    index=(vertices(k,1:3))';
    % riga sotto: opzione per funzione nota ai vertici
    uu_tmp=[uht(t_size,index), uht(t_size,index(1))];
    % riga sotto: opzione per funzione costante a tratti
    % p_tmp = ph(k)*ones(length(index)+1,1); 
    vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
    fill3(vert_temp(:,1),vert_temp(:,2),uu_tmp,uu_tmp); 
end
view(3)
grid on
colorbar
hold off







    
    