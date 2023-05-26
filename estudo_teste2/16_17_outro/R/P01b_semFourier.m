clear all; close all; clc

% Constantes
m = 0.5; %Kg
I = 10^(-4); %Kgm^2
K = 5; %N/m
delta = 10^(-3); %Nm
epslon = 10 ^(-2); %N

% Vetores
h = 0.1;
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

j = 1;
for i = 600:850
    if ((z(i-1) < z(i)) && (z(i+1) < z(i)))
        ind(j) = i;
        maximos(j) = z(i);
        t_maximos(j) = t(i);
        j = j + 1;
    end
end

f_osci = length(maximos)/(t(850)-t(600));
% f_osci = 0.4800 Hz;

% Plots
figure(1)
plot(t(600:850),z(600:850), t_maximos, maximos, 'o');
xlabel('t (s)')
ylabel('z (m)')
