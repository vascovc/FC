clear all; close all; clc

% Assim, temos E(n = 1) = 19.4988 (da al�nea a), E(n = 2) = 55.5165, E(n = 3) = 87.5036, E(n = 4) = 115.5000

En =  [19.4988, 55.5165, 87.5036, 115.500];
ni = 1:4;

p = polyfit(ni,En,2) % a = -2.0053 b = 42.0257 c = -20.5196

n = 1:15;
En_ply = p(1)*n.^2 + p(2)*n + p(3);

figure(1)
plot(n,En_ply)
xlabel('n')
ylabel('En')

En_max = max(En_ply)
ultimo_n = find(En_ply == En_max); %ultimo_n = 10 Ent�o existem 10 estados ligados e E( n = 10 = En_max = 199.2049

% Constantes
b = 10;

% Pr�- aloca��es
h = 0.01;
x = 0:h:b;
Nx = length(x);

psi = zeros(1,Nx);
psi_prog = zeros(1,Nx);
psi_reg = zeros(1,Nx);

g = zeros(1,Nx);
V = 200*(1-exp(-2*(x-1))).^2;

i_match = round(7*Nx/11); 

% Vetores para o shooting
Ne = 1000;
E = zeros(1,Ne);
ddpsi = zeros(1,Ne);

% guesses
E(1) = En_max-0.5;
E(2) = En_max+0.5;

% Condi��es fronteira
psi(1) = 0;
psi(2) = h*10^(-3); % previs�o necess�ria pela f�rmula do m�todo de Numerov
psi(Nx-1) = h*10^(-3); % previs�o necess�ria pela f�rmula do m�todo de Numerov
psi(Nx) = 0;

psi_prog(1) = 0;
psi_prog(2) = h*10^(-3); % previs�o necess�ria pela f�rmula do m�todo de Numerov
psi_prog(Nx-1) = h*10^(-3); % previs�o necess�ria pela f�rmula do m�todo de Numerov
psi_prog(Nx) = 0;

psi_reg(1) = 0;
psi_reg(2) = h*10^(-3); % previs�o necess�ria pela f�rmula do m�todo de Numerov
psi_reg(Nx-1) = h*10^(-3); % previs�o necess�ria pela f�rmula do m�todo de Numerov
psi_reg(Nx) = 0;


% Goal do shooting
B = 0; % ddpsi = 0

tol = 10^(-3);

k = 1;
comparator = 1;
while abs(B - comparator) > tol
    
    g(1:Nx) = 2*(E(k)-V(1:Nx));
    
    % M. Numerov
    for i = 2:Nx-1
        psi_prog(i+1) = (1+h^2*g(i+1)/12)^(-1)*(-(1+h^2*g(i-1)/12)*psi_prog(i-1) + 2*(1-5*h^2*g(i)/12)*psi_prog(i));
    end
    
    psi_prog_match = psi_prog(i_match);
    dpsi_prog = ((-25/12)*psi_prog(i_match) + 4*psi_prog(i_match+1) - 3*psi_prog(i_match+2) + (4/3)*psi_prog(i_match+3)-(1/4)*psi_prog(i_match+4))/h;
    
    for i = Nx-1:-1:2
        psi_reg(i-1) = (1+h^2*g(i-1)/12)^(-1)*(2*(1-5*h^2*g(i)/12)*psi_reg(i)-(1+h^2*g(i+1)/12)*psi_reg(i+1));
    end
    
    psi_reg_match = psi_reg(i_match);
    psi_reg = psi_reg*psi_prog_match/psi_reg_match;
    
    dpsi_reg = ((25/12)*psi_reg(i_match) - 4*psi_reg(i_match-1) + 3*psi_reg(i_match-2) - (4/3)*psi_reg(i_match-3)+(1/4)*psi_reg(i_match-4))/h;
    psi(1:i_match) = psi_prog(1:i_match);
    psi(i_match:Nx) = psi_reg(i_match:Nx); 
  
    comparator = (dpsi_prog/psi_prog_match - dpsi_reg/psi_reg_match)/(dpsi_prog/psi_prog_match + dpsi_reg/psi_reg_match);
    ddpsi(k) = comparator;
  
    
    if k > 1
        if ddpsi(k)-ddpsi(k-1) == 0
            break
        else
            m = (ddpsi(k)-ddpsi(k-1))/(E(k)-E(k-1));
            E(k+1) = E(k) + (B-ddpsi(k))/m;
        end
    end
    k = k + 1;
end

% Normaliza��o
C = trapz(x,abs(psi).^2);
psi_norm = psi/sqrt(C);

figure(2)
plot(x,psi_norm)
xlabel('x')
ylabel('psi(norm)')

E(k)

% E(n = 10) = 199.5026