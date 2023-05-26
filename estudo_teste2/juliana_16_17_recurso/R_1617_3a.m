clear all, close all, clc

% Pré-alocações 
D = 1;
L = 2; 

t0 = 0; %s
tf = 5; %s
x0 = 0; %s
xf = L;

dt = 0.01; %s
dx = 0.1; %cm
t = t0:dt:tf;
x = x0:dx:xf; 

Nt = length(t);
Nx = length(x);

rho = zeros(Nx, Nt);

% Condições iniciais
rho(:,1) = 2*sin(pi*x/L);

% Condições fronteira
rho(1,:) = 0;
rho(Nx,:) = 0;

% Método de Crank-Nicolson

ETA = dt*D/dx^2;

A = (2/ETA+2)*eye(Nx-2);
for  i= 1:Nx-2
    if (i > 1)
        A(i,i-1) = -1;
    end
    if (i < Nx-2)
        A(i,i+1) = -1;
    end
end

b = zeros(Nx-2,1);

for j = 1: Nt -1
    b = rho(1:Nx-2,j)+(2/ETA-2)*rho(2:Nx-1,j)+rho(3:Nx,j);
    b(1) = b(1) + rho(1,j+1);
    b(Nx-2) = b(Nx-2) + rho(Nx,j+1);
    rho(2:Nx-1,j+1) = linsolve(A,b);
    
%     plot(rho(:,j))
%     pause(0.1)
end


figure(1)
mesh(t,x,rho)
figure(2)
contourf(x,t,rho')
colorbar
