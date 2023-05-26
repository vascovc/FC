clear all; close all; clc

% Constantes físicas
a = 0.0025; %m
l = 0.0027; %m
alfa = pi/4;

% Condições iniciais
S0 = 0;
X0 = 0;

% Pré-alocação dos vetores a usar na integração
theta0 = 0;
thetamax = pi/2 - alfa;
h = 0.0018;
theta = theta0:h:thetamax;

S = zeros(1,length(theta));
X = zeros(1,length(theta));
S(1) = S0;
X(1) = X0;

% Método de Euler
for i = 1:length(theta)-1
    Q = (a/l^2)*(X(i)+h/a)-sin(theta(i))/(a*S(i));
    W = (a/l)^2*(X(i)+h/a)-sin(theta(i))/Q;
    
    if theta(i) == 0
        dS = 2*l^2/(a*h);
        dX = 0;
    else
        dS = cos(theta(i))/W;
        dX = sin(theta(i))/W;
    end
    
    S(i+1) = S(i) + dS*h;
    X(i+1) = X(i) + dX*h;
end
