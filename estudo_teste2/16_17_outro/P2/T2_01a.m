clear all, close all, clc

% Constantes do sistema
D = 0.02;
m = 5;
L = 1;

% Pré-alocações
x0 = 0;
xf = L;
t0 = 0;
tf = 3;

dx = 0.01;
dt = 0.001;

x = x0:dx:xf;
t = t0:dt:tf;

Nx = length(x);
Nt = length(t);

u = zeros(Nx,Nt);

% Condições iniciais
u(:,1) = 0.5;
u((x<0.6),1) = 0;
u((x>0.8),1) = 0;

% Condições fronteira
u(1,:) = 0;
u(Nx,:) = 0;

% Método de Euler
for j = 1:Nt-1
    for i = 2:Nx-1
        u(i,j+1) = u(i,j) + dt*(D*(u(i-1,j) - 2*u(i,j) + u(i+1,j))/dx^2 + m*u(i,j)*(1-u(i,j)));
    end
    
%     plot(x,u(:,j))
%     pause(0.1)
end

figure(1)
mesh(t,x,u)
xlabel('t')
ylabel('x')
zlabel('u')
figure(2)
contourf(x,t,u')
colorbar

