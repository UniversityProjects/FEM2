clear all
close all
clc

% h_max Array Definition
h = [0.2, 0.15, 0.1, 0.075, 0.05];


for i=1:5
    
    %%% Meshes Building
    disp(['Mesh: ', num2str(i), '   - h_max = ', num2str(h(i))]);
    [xv,yv,vertices,edges,endpoints,boundary,boundedges] = meshgen(h(i));
    nver = length(xv);
    nedge = length(endpoints);
    nele = length(edges);
    n_p2p0(i) = 2*nver + 2*nedge + nele;
    n_mini(i) = 3*nver;
    
    
    %%% Error Computation
    disp(['Errors: i = ', num2str(i)]);
    [errL2_P2P0(i), errH1_P2P0(i), uh, ph] = P2P0(xv,yv,vertices,edges,endpoints,boundary,boundedges);
    [errL2_MINI(i), errH1_MINI(i), uh, ph] = MINI(xv,yv,vertices,edges,endpoints,boundary,boundedges);
    
end


figure();
loglog (n_p2p0, errL2_P2P0, '-*b');
hold on;
loglog (n_mini, errL2_MINI, '-*r');
title ('ErrL2');
legend ('P2P0','MINI');

figure();
loglog (n_p2p0, errH1_P2P0, '-*b');
hold on;
loglog (n_mini, errH1_MINI, '-*r');
title ('ErrH1');
legend ('P2P0','MINI');

