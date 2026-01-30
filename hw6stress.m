
clc;

sigma = [40 10 2.5; 10 20 -5; 2.5 -5 -10];

[eig_vector eig_val] = eig (sigma);

% sigma_1 = 44.145 kPa
%sigma_2 = 17.009 kPa
%sigma_3 = -11.154 kPa

%% problem 2
n = log(0.5)/(log(0.4/0.6));

A = (30/185)/(0.4^1.71);


