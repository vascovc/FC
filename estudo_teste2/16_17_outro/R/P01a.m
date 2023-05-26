clear all; close all; clc

% Constantes
m = 0.5; %Kg
I = 10^(-4); %Kgm^2
K = 5; %N/m
delta = 10^(-3); %Nm
epslon = 10 ^(-2); %N

% Vetores
h = 0.01;
t0 = 0;
tf = 100;
t = t0:h:tf;
Nt = length(t); 

z = zeros(1,Nt);
theta = zeros(1,Nt);
v = zeros(1,Nt);
w = zeros(1,Nt);
Z = zeros(4,Nt);

% Condições iniciais
z(1) = 0.1; %m
v(1) = 0; 
theta(1) = 0;
w(1) = 0;

Z(1,:) = z;
Z(2,:) = v;
Z(3,:) = theta;
Z(4,:) = w;

% Matriz A
A =zeros(4,4);
A(1,1) = (h*K)/(2*m);
A(1,2) = 1;
A(1,3) = (h*epslon)/(4*m);
A(2,1) = 1;
A(2,2) = -h/2;
A(3,1) = (h*epslon)/(2*I);
A(3,3) = (h*delta)/(2*I);
A(3,4) = 1;
A(4,3) = 1;
A(4,4) = -h/2;

% Método de Crank-Nicolson
for i = 1:Nt-1
    Z(:,i+1) = linsolve(A,[Z(2,i) - ((h*K)/(2*m))*Z(1,i)- ((h*epslon)/(4*m))*Z(3,i); Z(1,i) + (h/2)*Z(2,i); Z(4,i) - ((h*epslon)/(4*I))*Z(1,i) - ((h*delta)/(2*I))*Z(3,i); Z(3,i) + (h/2)*Z(4,i)]);
end

z = Z(1,:);
v = Z(2,:);
theta = Z(3,:);
w = Z(4,:);

% Plots
figure(1)
plot(t,z);
xlabel('t (s)')
ylabel('z (m)')
figure(2)
plot(t,theta);
xlabel('t (s)')
ylabel('theta (rad)')

figure(3)
plot(t,v);
xlabel('t (s)')
ylabel('v (m/s)')
figure(4)
plot(t,w);
xlabel('t (s)')
ylabel('w (rad/s)')
% Observando as figuras 3 r 4 observa-se que quando w é máximo, v é mínimo
% e quando v é máximo, w é mínimo, o que vai de encontro ao facto do
% pêndulo alternar entre oscilações puramente longitudinais e oscilações
% puramente torsionais

figure(5)
plot(t,z,t,v)
% Observando o plot da figura 5 verifica-se que quando z = 0, v = 0, o que
% significa que quando Em = 0, o pendulo está em z = 0 e tem v = 0.

figure(6)
plot(t,theta,t,w)
% Observando o plot da figura 6 verifica-se que quando theta = 0, w = 0

figure(7)
plot(t,v/max(v),t,w/max(w))
% A figura 7 sustenta o que se verifica pela observação das figuras 3 e 4