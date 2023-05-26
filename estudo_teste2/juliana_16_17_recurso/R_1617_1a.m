clear all; close all; clc

%ODE45 constants
reltol = 4*10^(-13);
abstol_1 = 1*10^(-13);
abstol_2 = 1*10^(-13);
abstol_3 = 1*10^(-13);

options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2 abstol_3]);

% Constantes físicas
B = [0.3 0.4 0.6 0.72];

% Condições iniciais
t0 = 0;
tf = 50;
x0 = 0.1;
v0 = 0.1;
a0 = 0;


for k = 1:length(B)
    [t,sol] = ode45(@f,[t0 tf],[x0 v0 a0],options,B(k));
    
    x = sol(:,1);
    v = sol(:,1);
    a = sol(:,1);
    
    figure(1)
    plot(t,x)
    figure(2)
    plot(t,v)
    figure(3)
    plot(x,v)
end
