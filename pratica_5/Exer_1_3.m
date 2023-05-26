clear all
close all
%solucao completa
u = 10^-3; L = 1; T = 10^3;
n = [1 2 3];
wn = n*pi/L*sqrt(T/u);
x0 = 0;
h = 0.01;
x = x0:h:L;
N = length(x);
v = nan(1,N);
y = nan(1,N);
y(1) = 0;
v(1) = 2E-2;

guess = [2000 3500];
result = nan(1,2);

for i = 1:1:2
    for k = 1:N-1
        v(k+1) = v(k) + -guess(i)^2*u/T*y(k) * h;
        y(k+1) = y(k) + v(k+1) * h;
    end
    result(i) = y(N);
end

m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
b = result(2)-m*guess(2);
B = 0;%neste caso, quer-se 0
guess(3) = guess(2) + (B-result(2))/m;

%para obter a melhor estimativa
while abs( guess(2)-guess(1) ) > 1E-6
    %para manter o ciclo limpo
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    
    %Euler-Cromer
    for k = 1:N-1
        v(k+1) = v(k) + -guess(2)^2*u/T*y(k) * h;
        y(k+1) = y(k) + v(k+1) * h;
    end
    result(2) = y(N);
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    %b = result(2)-m*guess(2); % nao e necessario
    guess(3) = guess(2) + (B-result(2))/m;
    %seria de recomendar usar um valor limite para ver a covergencia
end

disp(['solucao obtida    ',num2str(guess(3))])
disp(['solucao analitica ',num2str(wn(1))])