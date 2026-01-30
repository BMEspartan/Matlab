% close all
clear variables
clc
% problem 1
k0 = 700; % pN/um
k1 = 2000/1.6; % pN/um
eta0 = 150; % pN*s/um
tau = eta0*(k0+k1)/(k0*k1); % s
F0 = 2000; % pN

t = linspace(0,2,100); % s
x = @(t) F0/k1*(1-k0/(k0+k1).*exp(-t/tau)); % um

plot(t,x(t)); % plot
ylim([0,2]);% duration of first rise 2 seconds
xlabel('Time(s)', 'FontSize',12);
ylabel('Displacement x(t) (um)', 'FontSize',12); 
%% 
% problem 2 
k0 = 700; % pN/um
k1 = 2000/1.5; % pN/um
eta0 = 150; % pN*s/um
tau = eta0*(k0+k1)/(k0*k1); % s
F0 = 2000; % pN

t = linspace(10,12,100); % s
x = @(t) F0/k1*(1-k0/(k0+k1).*exp(-(t-10)/tau)); % um

plot(t,x(t));
ylim([0,2]);
xlabel('Time(s)');
ylabel('Displacement x(t) (um)');