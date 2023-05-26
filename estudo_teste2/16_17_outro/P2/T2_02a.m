clear all, close all, clc

% Constantes do sistema
L = 2.00;
c = 1.50;
a = 1.00 * 10^(-3);
C = 1;

% Pré-alocações
x0 = 0;
xf = L;
t0 = 0;
tf = 7;

dx = 0.1;
dt = C*dx/c;

x = x0:dx:xf;
t = t0:dt:tf;

Nx = length(x);
Nt = length(t);

y = zeros(Nx,Nt);
C = c*dt/dx;

% Condições iniciais
y(:,1) = 0;
pd = a*x.*(L-x);

% Condições fronteira
y(1,:) = 0;
y(Nx,:) = 0;

% Resolução numérica
for n = 1: Nt-1
    for i = 2:Nx-1
        if n ==1
            y(i,2) = (C^2/2)*(y(i-1,1)-2*y(i,1)+y(i+1,1)) + y(i,1) + dt*pd(i);
        else
            y(i,n+1) = C^2*(y(i-1,n)-2*y(i,n)+y(i+1,n)) + 2*y(i,n) - y(i,n-1);
        end
    end
    
%     plot(x,y(:,n))
%     axis([x0, xf, min(y(:)), max(y(:))])
%     pause(0.2)
    
end

figure(1)
mesh(t,x,y)
xlabel('t')
ylabel('x')
zlabel('y')
figure(2)
contourf(x,t,y')
colorbar

