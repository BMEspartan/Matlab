%% Problem 1
k0=750; % pN/um
k1=2000/1.6; % pN/um
eta0=187.5; % pN*s/um


tau=eta0*(k0+k1)/(k0*k1);

F0=2000;% pN
t=0:.001:2;
x=F0*(1/k1)*(1-(k0/(k0+k1))*exp(-t/tau));
figure(1)
plot(t,x)
ylabel("Disp (um)")
xlabel("TIme (s)")

%% Problem 2
% k0=4000/1.5;
% k1=2000/1.5;
% eta0=533;
k0=750;
k1=2000/1.5;% just changed k1
eta0=187.5;

tau=eta0*(k0+k1)/(k0*k1);

F0=2000;
t=0:.001:2;
x=F0*(1/k1)*(1-(k0/(k0+k1))*exp(-t/tau));
figure(2)
plot(t,x)
ylabel("Disp (um)")
xlabel("TIme (s)")