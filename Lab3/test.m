clc
clear all

% Eulero all'indietro
thetaEI = 1;

% Eulero all'indietro
thetaEE = 0;

% CN
thetaCN = 0.5;

% i
t_step = [5, 10, 20, 40, 80];

% j
h_max = [0.2, 0.1, 0.05];


for i=1:5
    for j=1:3
        disp(['i = ', num2str(i)]);
        disp(['j = ', num2str(j)]);
        
        errHei = parabolic(h_max(j),t_step(i),thetaEI);
        ErrH1ei(i,j) = errHei;
        
        errHee = parabolic(h_max(j),t_step(i),thetaEE);
        ErrH1ee(i,j) = errHee;
        
        errHcn = parabolic(h_max(j),t_step(i),thetaCN);
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

