clc
clear all
close all

% out = 'yes';
out = 'no';

% plot = 'yes';
plot = 'no';



% h_max Array Definition
h = [0.800, 0.400, 0.200, 0.150, 0.100, 0.075, 0.050];
n = length(h);

% Variables Definition
% uh = zeros(n);
errL2 = zeros(1,n);
errH1 = zeros(1,n);

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
   [uh, errL2(i), errH1(i)] = main_P2_Vect (...
                                                xv, yv, ...
                                                vertices, ...
                                                edges, ...
                                                endpoints, ...
                                                boundary, ...
                                                boundedges ...
                                                );

end

clc

disp(['h = ', num2str(h)]);
disp(' ');

% Plot Loop
errL2
errH1
    
% Figure 1
figure(10);
hold on;
loglog (h, errL2, '-*r', h, h.^3,'-b');
grid on;
legend ('ErrL2', 'h^3', 'location', 'northeastoutside');
title ('ErrL2');
hold off;
  
% Figure 2
figure(11);
hold on;
loglog (h, errH1, '-*r', h, h.^2,'-b');
grid on;
legend ('ErrH1', 'h^2', 'location', 'northeastoutside');
title ('ErrH1');
hold off;