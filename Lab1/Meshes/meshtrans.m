function [xv,yv,vertices,edges,endpoints,boundary,boundedges] = meshtrans(P,E,T)
% --------------------------------------------------
% Prende le matrici P,E,T generate con PDETOOL ...
% ... e le trasforma nel formato nostro.
% Vedi meshgen.m per informazioni sulla struttura dati
%---------------------------------------------------

P=P';
T=T';
m=size(T,1);
T=[3*ones(m,1),T(:,1:3)];

mesh = struct ('V',P,'P',T);
mesh.E = edgefinder(T);
mesh.B = bvfinder(mesh.E); 

%-----------
% Genero mesh.PE

M=sparse(zeros(size(mesh.V,1)));
for k=1:size(mesh.E,1)    % build auxiliary matrix M
    M(mesh.E(k,1),mesh.E(k,2))=k;
end    
M=M+M';
for k=1:size(mesh.P,1)
    for j=1:mesh.P(k,1)-1
        mesh.PE(k,j)=M(mesh.P(k,j+1),mesh.P(k,j+2));
    end
    mesh.PE(k,mesh.P(k,1))=M(mesh.P(k,mesh.P(k,1)+1),mesh.P(k,2));
end

%-----------
% Infine converto lo struct "mesh" nel formato nostro:

xv = mesh.V(:,1);
yv = mesh.V(:,2);
vertices = mesh.P(:,2:end);
edges = [mesh.PE(:,2:3),mesh.PE(:,1)];  % ri-ordinamento colonne 
edges = full(edges);
boundary = mesh.B;
endpoints = mesh.E(:,1:2);

%-----------
% Genero boundedges a parte,sfruttando mesh.E (vedi VEM codes)

boundedges=[];
for i=1:size(mesh.E,1)
    if mesh.E(i,4)==0
        boundedges=[boundedges,i];
    end    
end


% numele         % opzione per stampare numeri di elemento sulla mesh  
% numver         % opzione per stampare numeri di vertice sulla mesh
% numedge        % opzione per stampare numeri di lato sulla mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------
function [E]=edgefinder(T)
%-----------------------------------------------------------
% given a mesh defined by input matrix T it
% finds the edges and generates an edge matrix E organized as follows.
% On each row there appear the vertex number 1 number, vertex number 2 number, 
% ... and then left and right element numbers. 
% NOTE: matrix E has a zero in the last column instead of 
% ... an element number if this is a boundary edge.
% ---------
% Input matrix T has on each row
% the data for an element, i.e. the number of vertexes, then the vertex
% numbers ordered in anti-clockwise sense.
% ---------
m=size(T,1);
k=1;
% ----- 
% generate matrix A with repeated edges, generated element by element 
% -----
for i=1:m
    for j=1:(T(i,1)-1)
        A(k,1:2)=T(i,j+1:j+2);
        A(k,3)=i;
        k=k+1;
    end
    A(k,1)=T(i,T(i,1)+1);
    A(k,2)=T(i,2);
    A(k,3)=i;
    k=k+1; 
end       
% ----- 
%    generate matrix E working on matrix A
% -----
flag=0;
% N=size(A,1);
for k=1:size(A,1);
    flag=0;
    for kk=k+1:size(A,1)
        if (( A(k,1)==A(kk,2) ) && ( A(k,2)==A(kk,1) ))   
            E(k,:)=[A(k,1:3),A(kk,3)];
            A=[A(1:kk-1,:) ; A(kk+1:end,:)];    
            % N=N-1;    
            flag = 1;
            break
        end
    end
    if k > size(A,1)  
       break
    end
    if flag < 0.5
        E(k,:)=[A(k,1:3),0];
    end
end    


%-----------------------------------------------------------
function[V]=bvfinder(E);
%-----------------------------------------------------------
% given a mesh, this function needs in input the edge matrix generated by
% the function edgefinder. Then, it gives a vector which contains the 
% numbers of the nodes which are on the boundary
% -------
% Generate preliminary vector V with repetitions
% -------
k=1;
for i=1:size(E,1)
    if E(i,4)==0  % this determines boundary edges, see edgefinder function
        V(k*2-1) = E(i,1);
        V(k*2) = E(i,2);
        k=k+1;
    end
end
% ------
% Erase repetitions
% ------
for i=1:length(V)
    for j=i+1:length(V)
        if V(i)==V(j)
            V = [V(1:j-1),V(j+1:end)];
            break
        end        
    end
end
