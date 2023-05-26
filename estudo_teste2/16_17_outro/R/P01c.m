clear all; close all; clc

% Constantes
m = 0.5; %Kg
I = 10^(-4); %Kgm^2
K = 5; %N/m
delta = 10^(-3); %Nm
epslon = 10 ^(-2); %N

% Vetores
h = 0.001;
t0 = 0;
tf = 100;
t = t0:h:tf;
Nt = length(t); 

z = zeros(1,Nt);
theta = zeros(1,Nt);
v = zeros(1,Nt);
w = zeros(1,Nt);
Z = zeros(4,Nt);

Em = zeros(1, Nt);
Ep = zeros(1, Nt);

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

Em(1) = K*Z(1,1)^2/2 + m*Z(2,1)^2/2;
Ep(1) = K*Z(1,1)^2/2;

% Método de Crank-Nicolson
for i = 1:Nt-1
    Z(:,i+1) = linsolve(A,[Z(2,i) - ((h*K)/(2*m))*Z(1,i)- ((h*epslon)/(4*m))*Z(3,i); Z(1,i) + (h/2)*Z(2,i); Z(4,i) - ((h*epslon)/(4*I))*Z(1,i) - ((h*delta)/(2*I))*Z(3,i); Z(3,i) + (h/2)*Z(4,i)]);
    Em(i+1) = K*Z(1,i+1)^2/2 + m*Z(2,i+1)^2/2;
    %Ep(i+1) = K*Z(1,i+1)^2/2;
end

z = Z(1,:);
v = Z(2,:);
theta = Z(3,:);
w = Z(4,:);

% Plots
figure(1)
plot(t,Em);
xlabel('t (s)')
ylabel('Em (J)')

figure(2)
plot(t,Em/max(Em),t,v/max(v),t,w/max(w));

% Pela observação da figura 7 nota-se que o perfil da Energia Mecânica
% coincide com o perfil do envelope da velocidade. Assim, o intervalo de
% tempo entre oscilações puramente longitudinais e puramente torsionais
% ocorre entre v = 0 e v ? vmax, o que corresponde ao intervalo entre Em =
% 0 e Em = Emax

figure(3)
plot(t(27000:53000), Em(27000:53000)) % utilizando apenas estes pontos é mais fácil determinar máximos e mínimos
Em_min = min(Em(27000:53000));
Em_max = max(Em(27000:53000));
ind_max = find(Em == Em_max); 
ind_min = find(Em == Em_min); 
t_min = t(ind_min);
t_max = t(ind_max);
intervalo_tempo = t_max-t_min % Assim, o intervalo entre oscilações puramente longitudinais e oscilações puramente torsionais é de 11.3970s

