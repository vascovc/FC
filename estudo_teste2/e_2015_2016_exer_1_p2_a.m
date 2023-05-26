clear all
close all

beta = 18;
x0 = 0;
xf = 7;
h = 0.01;

x = x0:h:xf;
N = length(x);

T = nan(1,N);
Tv = nan(1,N);
T(1) = 10^-4;
Tv(1) = 10^-4;

B = 0;

guess = [1.9 2.0];
result = nan(1,2);

for i = 1:1:2
    for k = 1:N-1
        Tv(k+1) = Tv(k) + (guess(i)*Tv(k)-beta*exp(-1/T(k))*( 1+Tv(k)-guess(i)*T(k) ))* h;
        T(k+1) = T(k) + Tv(k) * h;
    end
    result(i) = Tv(N);
end

m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
guess(3) = guess(2) + (B-result(2))/m;

%para obter a melhor estimativa
while abs( guess(2)-guess(1) ) > 1E-4
    %para manter o ciclo limpo
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    
    %Euler
    for k = 1:N-1
        Tv(k+1) = Tv(k) + (guess(2)*Tv(k)-beta*exp(-1/T(k))*(1+Tv(k)-guess(2)*T(k) ))* h;
        T(k+1) = T(k) + Tv(k) * h;
    end
    result(2) = Tv(N);
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    %b = result(2)-m*guess(2); % nao e necessario
    guess(3) = guess(2) + (B-result(2))/m;
    %seria de recomendar usar um valor limite para ver a covergencia
end
subplot(1,2,1)
plot(x,T)
subplot(1,2,2)
plot(x,Tv)

