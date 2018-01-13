function [dat,helix] = helix_heteroscedastic;

%% DATA: NOISY HELIX

N = 10000;

a = 3;
b = 2;
t = linspace(0,4*pi,N);
x = a*cos(t);
y = a*sin(t);
z = b*t;

helix = [x; y; z];

ruido1 = [];
ruido2 = [];
sigma = linspace(0.1,1.8,length(x));
for i=1:length(x)
    ruido1 = [ruido1 sigma(i)*randn(3,1)]; 
    ruido2 = [ruido2 sigma(i)*randn(3,1)]; 
end

dat = helix + ruido1;
dat2 = helix + ruido2;
