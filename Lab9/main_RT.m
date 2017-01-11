function [uh,ph]=main_RT(xv,yv,vertices,edges,endpoints);
%-------------------------------
% PROGRAMMA PER GLI RT_0
% Prende in input i dati della mesh (vedi meshgen.m per descrizione)e
% ... rende in output flussi e pressioni. Inoltre fa varie stampe e 
% ... calcola vari errori (se e' nota soluzione esatta).
% -------------------------------
% Il dato di carico (sorgente) dato tramite la function carico.m
% La soluzione eventuale data tramite le function exsolF.m, exsolP.m
% Si considera caso con dirichlet omogenee sulla variabile scalare
% -------------------------------
flag=1;   % flag =1 calcola errori, flag=0 non calcola errori
mm = size(vertices,1);   % numero di elementi
ndof = size(endpoints,1);  % numero G. di Lib. Flussi(equivale al numero di lati)
S = sparse(ndof+mm,ndof+mm); % crea matrice vuota (formato sparso)  
F = zeros(ndof+mm,1);     % crea vettore vuoto (sara' per il carico)
elemdata = struct;   % serve a salvare dati utili per ogni elemento,
                     % ... che servono nel post-processing (calc errore)
% -------------
% Matrice di stiffness                     
% -------------        
for k=1:mm          % cicla sugli elementi     
    Xv = xv(vertices(k,1:3))';      % le x-coordinate dei 3 vertici
    Yv = yv(vertices(k,1:3))';      % le y-coordinate dei 3 vertici   
    [K,P,f,M,E] = RT_loc(Xv,Yv);    % matrici di stiffness locali    
    dofvec = edges(k,1:3);  % vettore degli indici di lato del elemento     
    % Sotto aggiusto cambi di segno dovuti alle direzioni delle normali
    % Nei GdL globali, le normali sono prese come segue: per ogni lato,
    % ... si prende il vettore che va dal vertice di numero inferiore
    % ... a quello di numero superiore (nel ordinamento globale dei 
    % ... vertici) e poi si ruota in senso orario di pi/2 (e si normalizza).
    r=ones(3,1);
    if vertices(k,3)-vertices(k,2) < 0
       r(1)=-1;
    end   
    if vertices(k,1)-vertices(k,3) < 0
       r(2)=-1;
    end    
    if vertices(k,2)-vertices(k,1) < 0
       r(3)=-1;
    end    
    H = diag(r); % matrice dei cambi di segno dovuti alle normali    
    S(dofvec,dofvec) = S(dofvec,dofvec) + H*K*H ; % parte flussi-flussi
    S(dofvec,ndof+k) = S(dofvec,ndof+k) + H*P; % flussi-press
    S(ndof+k,dofvec) = S(ndof+k,dofvec) + (H*P)'; % press-flussi
    F(ndof+k)=f;
    % salvo dati utili per dopo: 
    elemdata(k).MH = M*H;
    Areas(k) = E;
end 


%----------
sol = S\F;     % risoluzione sistema lineare
%----------
uh = sol(1:ndof);        % estraggo i flussi
ph = sol(ndof+1:end);    % estraggo le pressioni
clear sol;

% ---------------------------------------------------------
% POST-PROCESSING: calcola valori ai baricentri degli elementi
% ---------------------------------------------------------
for k=1:mm          % cicla sugli elementi
    MH = elemdata(k).MH;  
    dofvec = edges(k,1:3);    % vettore degli indici di lato del elemento     
    uloc = MH*uh(dofvec);     % calcola gli [a,b,c] locali sul elemento
    elemdata(k).uloc = uloc;  % salva [a,b,c] del elemento per uso post                        
    uhbar(k,:) = uloc(1:2);   % il valore nel baricentro sono [a,b]           
end 


if flag==1  % se flag=0 non calcola gli errori   
    
% INTERPOLAZIONE AI BARICENTRI DI FLUSSI E PRESSIONI ESATTE:
uI=zeros(mm,2);   % flussi esatti interpolati, una colonna per componente
pI=zeros(mm,1);   % pressioni esatte interpolate 
for k=1:mm          % cicla sugli elementi     
    Xv = xv(vertices(k,1:3));      % le x-coordinate dei 3 vertici
    Yv = yv(vertices(k,1:3));      % le y-coordinate dei 3 vertici 
    b = ((sum([Xv,Yv]))')/3;       % baricentro del elemento
    [fex]=exsolF(b);
    [pex]=exsolP(b);
    uI(k,:) = fex;
    pI(k) = pex;
end  
 
%---------------------------------------------------------
%  ERRORI VARI (TUTTI RELATIVI)
% --------------------------------------------------------

% ERRORE IN NORMA L2 GREZZO (APPROX INTEGRALE USANDO BARICENTRI)
% Areas: vettore riga con le aree di tutti gli elementi (ordinate)
ss = sum((uhbar'-uI').^2)';
disp('errore L2 flussi solo ai baricentri')
errfluxb = sqrt( (Areas*ss)/(Areas*( (sum((uI').^2))' ) ) )
ss = sum((ph'-pI').^2)';
disp('errore L2 pressioni solo ai baricentri')
errpressb = sqrt( (Areas*ss)/(Areas*( (sum((pI').^2))' ) ) )


% ERRORE IN NORMA L2 PIU FINE (FORMULA INT. CON PUNTI MEDI)
errpress2 = 0;  % errore pressioni assoluto al quadrato
press2 = 0;  % norma L2 delle pressioni
errflux2  = 0;  % errore flussi assoluto al quadrato
flux2  = 0;  % norma L2 dei flussi
for k=1:mm          % cicla sugli elementi   
   uloc = elemdata(k).uloc;   % carico le [a,b,c] locali salvate prima      
   % MH = elemdata(k).MH;     
   % dofvec = edges(k,1:3);  % vettore degli indici di lato del elemento     
   % uloc = MH*uh(dofvec);    % calcola gli [a,b,c] locali sul elemento            
   Xv = xv(vertices(k,1:3))';      % le x-coordinate dei 3 vertici
   Yv = yv(vertices(k,1:3))';      % le y-coordinate dei 3 vertici   
   X = [Xv;Yv];
   P = [X(:,2)+X(:,3),X(:,3)+X(:,1),X(:,1)+X(:,2)]/2; % matrice dei punti medi     
   errpress2 = errpress2 + Areas(k)*( (exsolP(P(:,1))-ph(k))^2 + ...
               (exsolP(P(:,2))-ph(k))^2 + ...
               (exsolP(P(:,3))-ph(k))^2 )/3;
   press2 = press2 + Areas(k)*( exsolP(P(:,1))^2 + ...
               exsolP(P(:,2))^2 + exsolP(P(:,3))^2 )/3;
   b = ((sum(X'))')/3;     % baricentro del elemento
   % calcolo uhloc, che sono i flussi discreti nei tre punti medi,
   % ... organizzati come matrice 2x3 (colonna 1 ha x e y nel primo punto medio)
   uh_ptimed = [uloc(1:2),uloc(1:2),uloc(1:2)] + uloc(3)*(P - [b,b,b]); 
   errflux2 = errflux2 + Areas(k)*( ...
                norm((exsolF(P(:,1))-uh_ptimed(:,1)))^2 + ...
                norm((exsolF(P(:,2))-uh_ptimed(:,2)))^2 + ...
                norm((exsolF(P(:,3))-uh_ptimed(:,3)))^2 )/3;  
   flux2 = flux2 + Areas(k)*( ...
                norm(exsolF(P(:,1)))^2 + ...
                norm(exsolF(P(:,2)))^2 + ...
                norm(exsolF(P(:,3)))^2 )/3;           
end 
disp('errore L2 pressioni')
errpress = sqrt(errpress2/press2)
disp('errore L2 flussi')
errflux = sqrt(errflux2/flux2)
end       % end del "if" con flag


%---------------------------------------------------------
%  PLOT DI FIGURE VARIE
% --------------------------------------------------------

% PLOT DELLE PRESSIONI
figure(1)
for k=1:size(vertices,1)
    hold on;
    index=(vertices(k,1:3))';
    % riga sotto: opzione per funzione nota ai vertici
    % uu_tmp=[uu(index);uu(index(1))];
    % riga sotto: opzione per funzione costante a tratti
    p_tmp = ph(k)*ones(length(index)+1,1); 
    vert_temp=[xv(index),yv(index); xv(index(1)),yv(index(1))];
    fill3(vert_temp(:,1),vert_temp(:,2),p_tmp,p_tmp); 
end
view(3)
grid on
colorbar
hold off

% PLOT DEI FLUSSI
% (plotta i flussi come vettori nei baricentri degli elementi)
TP = []; 
for k=1:mm          % cicla sugli elementi   
   uloc = elemdata(k).uloc;   % carico le [a,b,c] locali salvate prima  
   Xv = xv(vertices(k,1:3))';      % le x-coordinate dei 3 vertici
   Yv = yv(vertices(k,1:3))';      % le y-coordinate dei 3 vertici   
   X = [Xv;Yv];
   b = ((sum(X'))')/3;   % baricentro del elemento
   TP = [TP; b(1) , b(2), uloc(1), uloc(2)]; %
end
figure(2)
quiver(TP(:,1),TP(:,2),TP(:,3),TP(:,4));
hold off;
   
   

