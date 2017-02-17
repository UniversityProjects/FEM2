clc
clear all
close all

% out = 'yes';
out = 'no';

% plot = 'yes';
plot = 'no';

% Lambda Array Definition
% lambda = [0.0001, 0.001, 1, 100, 10000, 1000000];
lambda = [0.0001, 1, 1000, 1000000, 10000000000];
nl = length(lambda);

% h_max Array Definition
h = [0.4, 0.2, 0.15, 0.1, 0.075, 0.05];
n = length(h);

% Variables Definition
% uh = zeros(n);
errL2 = zeros(n,nl);
errH1 = zeros(n,nl);

for j=1:nl
    
    disp(['lambda = ', num2str(lambda(j))]);
    
    carico(lambda(j));
    
for i=1:n
    
    %%% Meshes Building
    disp(['    Mesh: ', num2str(i), '   - h_max = ', num2str(h(i))]);
    meshname =  num2str(h(i));
    [xv,yv,vertices,edges,endpoints,boundary,boundedges] = meshgen(h(i));
    
    nver = length(xv);
    nedge = length(endpoints);
    nele = length(edges);
   % n_p2(i) = 2*nver + 2*nedge + nele;
   % meshplot(xv,yv,endpoints,meshname,i);
    
    
    %%% Error Computation
    % [uh(i), errL2(i), errH1(i)] = main_P2_Vect (...
    [errL2(i,j), errH1(i,j)] = P2_Vect (...
                                         xv, yv, ...
                                         vertices, ...
                                         edges, ...
                                         endpoints, ...
                                         boundary, ...
                                         boundedges, ...
                                         lambda(j), ...
                                         out, plot ...
                                         );
                                     
    disp(['    errL2 = ', num2str(errL2(i,j))]);
    disp(['    errH1 = ', num2str(errH1(i,j))]);
                                   
    
end
end

clc

disp(['h = ', num2str(h)]);
disp(' ');

% Plot Loop
 for j=1:nl
     
  disp(['lambda = ', num2str(lambda(j))]); 
  disp(['     errL2 = ', num2str(errL2(:,j)')]);
  disp(['     errH1 = ', num2str(errH1(:,j)')]); 
    
% Figure 1
  figure(1);
  hold on;
 
    subplot(3, 2, j);
    loglog (h, errL2(:,j), '-*r', h, h.^3,'-b');
    grid on;
    % legend ('ErrL2', 'h^3', 'location', 'northeastoutside');
    title (['ErrL2 (L = ', num2str(lambda(j)), ')']);
  
  saveas (1, 'convL2.png');
  
% 
%  Figure 2
  figure(2);
 
  hold on;
  
      subplot(3, 2, j), loglog (h, errH1(:,j), '-*r', h, h.^2,'-b');
      grid on;
      % legend ('ErrH1', 'h^2', 'location', 'northeastoutside');
      title (['ErrH1 (L = ', num2str(lambda(j)), ')']);

  saveas (2, 'convH1.png'); 
  
     
  end

  
 % LATEX ERROR TABLE
 latex = 'yes';
 
if (strcmp(latex,'yes'))
 for j=1:nl
 
 sL2 = 'ErrL2	&	$';
 sH1 = 'ErrH1	&	$';   
     
  disp(['lambda = ', num2str(lambda(j))]); 
  
  for i=1:n
      sL2 = [sL2, num2str(errL2(i,j))];
      sH1 = [sH1, num2str(errH1(i,j))];
      
      if (i==n)
         sL2 = [sL2, '$ \\ \hline'];
         sH1 = [sH1, '$ \\ \hline']; 
      else
          sL2 = [sL2, '$ & $'];
          sH1 = [sH1, '$ & $'];         
      end
  end
  
  disp(sL2);
  disp(sH1);
  disp('');
  
 end 
  
end % End If
  

  

