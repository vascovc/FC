clear all, close all, clc

% Pré-alocações 
D = 1;
L = 2; 

t0 = 0; %s
tf = 50; %s
x0 = 0; %s
xf = L;

dt = 0.1; %s
dx = 0.1; %cm
t = t0:dt:tf;
x = x0:dx:xf; 

Nt = length(t);
Nx = length(x);

T = zeros(Nx, Nt);

% Condições iniciais
T(:,1) = 2*sin(pi*x/L);

% Condições fronteira
T(1,:) = 0;
T(Nx,:) = 0;

ETA = dt*D/dx^2;

% Método de Crank-Nicolson
A = (2/ETA + 2)*eye(Nx-2);
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
    for i = 2:Nx-1
        b(i-1) = T(i,j) + (2/ETA-2)*T(i+1,j) + T(i+2,j);
        b(1) = b(1) + T(1,j+1);
        b(Nx-2) =  b(Nx-2) + T(Nx,j+1);
        T(i,j+1) = linsolve(A,b);
    end
end



figure(1)
mesh(t,x,T)
figure(2)
contourf(x,t,T')
colorbar