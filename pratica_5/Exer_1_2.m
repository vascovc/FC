clear all
close all
%primeira iteracao; ver se funcionava tudo
u = 10^-3; L = 1; T = 10^3;%dados
n = [1 2 3];
wn = n*pi/L*sqrt(T/u);% solucao analitica
x0 = 0;
h = 0.01;
x = x0:h:L;
N = length(x);
v = nan(1,N);%velocidade
y = nan(1,N);%posicao
y(1) = 0;
v(1) = 2E-2;%irrelevante

guess = [2000 3500];%guess inicial
result = nan(1,2);
for i=1:1:2
    for k = 1:N-1%metodo de Euler-Cromer
        v(k+1) = v(k) + -guess(i)^2*u/T*y(k) * h;
        y(k+1) = y(k) + v(k+1) * h;
    end
    result(i) = y(N);
end

m = ( result(2)-result(1) )/( guess(2)-guess(1) )
b = result(2)-m*guess(2)
guess(3) = guess(2) + (0-result(2))/m
