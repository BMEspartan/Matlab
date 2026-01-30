clc
clear variables
close all
v = @(s) 10.^(0.0871.*log10(s).^2-0.4823.*log10(s)+1.3167);
s = linspace(10^(-2),10^2,2000);
figure
loglog(s,v(s));
ylim([1,1000])
xlabel('Shear rate (sec^-1)', 'FontSize',12);
ylabel('Viscosity (cPoise) ', 'FontSize',12); 
figure
plot(s,v(s));
xlabel('Shear rate (sec^-1)', 'FontSize',12);
ylabel('Viscosity (cPoise) ', 'FontSize',12); 
%% 
% problem 2 c 

clc
clear variables
close all
p1 = 80*133.322; % mmHg to Pascal
miu = 3.5/100*0.1; % cP to kg/(m*s)
Q = 7.4/1000000; % mL/s to m^3/s
D = 6.967/1000; % mm to m
p2 = @(L) p1+128*miu*Q*L/(pi*D^4);
L = linspace(0,20/100,200);
plot(L*100,p2(L)/133.322) % 133.322
xlabel('Distance (cm)', 'FontSize',12);
ylabel('Pressure (mmHg)', 'FontSize',12); 
%% 
% problem 2 d
clc
clear variables
close all
miu = 3.5/100*0.1; % cP to kg/(m*s)
Q = 7.4/1000000; % mL/s to m^3/s
tau = @(D,miu) -32.*(miu/100*0.1).*Q./(pi.*D.^3);
D = linspace(1/1000,10/1000,200);
figure
plot(D*1000,tau(D,3.5)*10)
hold on
plot(D*1000,tau(D,4.5)*10)
hold on
plot(D*1000,tau(D,5.5)*10)
hold off
xlabel('Cycle-averaged diameter D', 'FontSize',12);
ylabel('Mean wall shear stress (dyn/cm^2)', 'FontSize',12); 

legend('3.5 cP', '4.5 cP', '5.5 cP',  'Location', 'best');