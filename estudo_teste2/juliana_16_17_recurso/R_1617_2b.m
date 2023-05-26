clear all; close all; clc

% Constantes físicas
a = 0.0025; %m
l = 0.0027; %m
alfa = pi/4;

% Condições iniciais
S0 = 0;
X0 = 0;

theta0 = 0;
thetamax = pi/2 - alfa;

% Pré-alocação dos vetores a usar na integração
N = 5000;
theta = zeros(1,N);
S = zeros(1,N);
X = zeros(1,N);
S(1) = S0;
X(1) = X0;
    
% Shooting goal
B = 1; %S(thetamax)

% Pré-alocação de vetores para o Shooting
N = 20; 
tol = 10^(-3);
h = zeros(1,N); %o que se quer descobrir
Smax = zeros(1,N); %o que se está a condicionar

% guess(1) e guess(2)
h(1) = 0.0018;
h(2) = 0.0019;

%Método de Shooting
comparator = 0; %para entrar no while
k = 1;
while abs(B - comparator) > tol
    n = floor(abs(theta0-thetamax)/h(k) + 1);
    theta(1:n) = theta0:h(k):thetamax;
    
    % Método de Euler
    for i = 1:n-1
        Q = (a/l^2)*(X(i)+h(k)/a)-sin(theta(i))/(a*S(i));
        W = (a/l)^2*(X(i)+h(k)/a)-sin(theta(i))/Q;
        
        if theta(i) == 0
            dS = 2*l^2/(a*h(k));
            dX = 0;
        else
            dS = cos(theta(i))/W;
            dX = sin(theta(i))/W;
        end
        
        S(i+1) = S(i) + dS*h(k);
        X(i+1) = X(i) + dX*h(k);
    end
    comparator = S(i+1);
    Smax(k) = S(i+1);
    
    if k > 1 
        declive = (Smax(k)-Smax(k-1))/(h(k)-h(k-1));
        h(k+1) = h(k) + (B-Smax(k))/declive;
    end
    
    k = k + 1;
end