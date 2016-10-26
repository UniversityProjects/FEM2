clc
clear all

% Eulero all'indietro
thetaEI = 1;

% Eulero all'indietro
thetaEE = 0;

% CN
thetaCN = 0.5;

% i
t_step = [5, 10, 20, 40, 80, 100];

% j
h_max = [0.5, 0.2, 0.1, 0.05];

% Matrices Definition
% ErrH1ei = zeros(6,4);
% ErrH1ee = zeros(6,4);
% ErrH1cn = zeros(6,4);

% Computation Loop
% Use evalc to suppress parabolic output
for i=1:6 % Time Loop
    for j=1:4 % Space Loop
        disp(['h = ', num2str(h_max(j)), '  ',...
              't = ', num2str(t_step(j)), '  ',... 
              'i = ', num2str(i), '  ',... 
              'j = ', num2str(j)]);
        
        [T, errHei] = evalc('parabolic(h_max(j),t_step(i),thetaEI);');
        ErrH1ei(i,j) = errHei;
        
        [T, errHee] = evalc('parabolic(h_max(j),t_step(i),thetaEE);');
        ErrH1ee(i,j) = errHee;
        
        [T, errHcn] = evalc('parabolic(h_max(j),t_step(i),thetaCN);');
        ErrH1cn(i,j) = errHcn;
    end
end

clc



clc

disp(['thetaEI = ', num2str(thetaEI)]);
ErrH1ei

disp(['thetaEE = ', num2str(thetaEE)]);
ErrH1ee

disp(['thetaCN = ', num2str(thetaCN)]);
ErrH1cn

