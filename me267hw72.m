k0= 750;
k1=2000/1.5;
eta0=1;

tau=eta0*(k0+k1)/ (k0*k1);

F0 = 200;
t= 0:.01:2;
x= F0 * ((1/k1) * (1- (k0/+k1)) *exp (-t/tau));
plot (t,x) 
ylabel ("Dispo (um) ")
xlabel ("Time (s) ")
