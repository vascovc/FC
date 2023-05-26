clear all, close all, clc

% Constantes do sistema
L = 2.00;
c = 1.50;
a = 1.00 * 10^(-3);
C = 1;

% Pr�-aloca��es
x0 = 0;
xf = L;
t0 = 0;
tf = 7.2;

dx_vec = [0.400 0.200 0.100 0.050 0.025];

NN = length(dx_vec);
erro_vec = zeros(1,NN);

dx_max = max(dx_vec);
dx_min = min(dx_vec);

Nx_vec = (xf-x0)./dx_vec + 1;

Nmin = (xf-x0)/dx_max + 1;
Nmax = (xf-x0)/dx_min + 1;

dt_vec = C*dx_vec./c;
dt_min = min(dt_vec);
Nt_vec = (tf-t0)./dt_vec +1;

x = x0:dx_min:xf; % guardar toda a mem�ria utilizada antes do for
t = t0:dt_min:tf; % guardar toda a mem�ria utilizada antes do for


Nx_max = length(x); % guardar toda a mem�ria utilizada antes do for
Nt_max = length(t); % guardar toda a mem�ria utilizada antes do for

y = zeros(Nx_max,Nt_max); % guardar toda a mem�ria utilizada antes do for

f = @(xi,ti) a*(2*c^2-xi*(L-xi))*sin(ti);


for idx = 1:NN
    dx = dx_vec(idx);
    dt = dt_vec(idx);
    Nx = Nx_vec(idx);
    Nt = Nt_vec(idx);
    
    x = x0:dx:xf;
    t = t0:dt:tf;
    
    % Solu��o anal�tica
    y_analitico = a*x.*(L-x)*sin(7.2);
    
    % Condi��es iniciais
    y(1:Nx,1) = 0;
    pd = a*x.*(L-x);
    
    % Condi��es fronteira
    y(1,1:Nt) = 0;
    y(Nx,1:Nt) = 0;

    % Resolu��o num�rica
    for n = 1: Nt-1
        for i = 2:Nx-1
            if n ==1
                y(i,2) = y(i,1) + (C^2/2)*(y(i-1,1)-2*y(i,1)+y(i+1,1)) + pd(i)*dt +(1/2)*dt^2*f(x(i),t(1));
            else
                y(i,n+1) = - y(i,n-1) + 2*y(i,n) + C^2*(y(i-1,n)-2*y(i,n)+y(i+1,n)) + dt^2*f(x(i),t(n));
            end
        end
    end
    
    erro_vec(idx) = max(abs(y_analitico(1:Nx)-y(1:Nx,Nt)'));
end

plot(log(dx_vec),log(erro_vec));
p = polyfit(log(dx_vec),log(erro_vec),1);

% O log do erro � proporcional ao log do dx, sendo o declive da reta m =
% 2.0052. Logo, deve existir uma rela��o: erro = const*dx^(2.0052). 