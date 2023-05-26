clear all; close all; clc

% Constantes
d = 0.010; %m
rho = 1.0 * 10^3; %Kg/m^3
mu = 1.0 * 10^(-3); %Pa*s

% Pré-alocação de vetores
t0 = 0; %s
tf = 100; %s
x0 = 0; %s
xf = d;

dt = 0.1; %s
dx = 0.0001; %m
t = t0:dt:tf;
x = x0:dx:xf; 

Nt = length(t);
Nx = length(x);

u = zeros(Nx, Nt);
G = 800*tanh(t/10);


% Condições iniciais
u(2:Nx-1,1) = 0;

% Condições fronteira
u(1,:) = 0;
u(Nx,:) = 0;

no = (mu*dt)/(rho*dx^2);

A = (2/no + 2)*eye(Nx-2);
b = zeros(Nx-2,1);

% Método de Crank-Nicolson
for  i= 1:Nx-2
    if (i > 1)
        A(i,i-1) = -1;
    end
    if (i < Nx-2)
        A(i,i+1) = -1;
    end
end

const = 2/no-2;

for j = 1: Nt -1
    b = u(3:Nx,j)+const*u(2:Nx-1,j)+u(1:Nx-2,j) + (dt/(no*rho))*(G(j) + G(j+1));
    b(1) = b(1) + u(1,j+1);
    b(Nx-2) = b(Nx-2) + u(Nx,j+1);
    u(2:Nx-1,j+1) = linsolve(A,b);
end

figure(1)
mesh(t,x,u)
xlabel('t')
ylabel('x')
figure(2)
contourf(x,t,u')
colorbar
