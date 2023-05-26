clear all
close all

a = 1;
h = 0.001;
x = -a:h:a;
N = length(x);

y = zeros(1,N);
y(1) = 0;
y(end) = 0;
y(2) = h/1000; %valor pequeno diferente de 0
n = 1;

g = zeros(1,N);
guess = [n^2*pi^2/(8*a^2)-1 n^2*pi^2/(8*a^2)+1];
result = nan(1,2);

for i = 1:1:2
    g = repmat(2*guess(i),N);
    for k = 2:N-1
        y(k+1) = (1+h^2/12*g(k+1))^-1*(-(1+h^2/12*g(k-1))*y(k-1)+2*(1-5*h^2/12*g(k))*y(k));
    end
    result(i) = y(N);
end
m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
b = result(2)-m*guess(2);
B = 0;%neste caso, quer-se 0 para a posicao final
guess(3) = guess(2) + (B-result(2))/m;

while abs( guess(2)-guess(1) ) > 1E-7
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    
    g = repmat(2*guess(2),N);
    for k = 2:N-1
        y(k+1) = (1+h^2/12*g(k+1))^-1*(-(1+h^2/12*g(k-1))*y(k-1)+2*(1-5*h^2/12*g(k))*y(k));
    end
    result(2) = y(N);
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    %b = result(2)-m*guess(2); % nao e necessario
    guess(3) = guess(2) + (B-result(2))/m;
    %seria de recomendar usar um valor limite para ver a covergencia
end

disp(['valor pratico obtido ',num2str(guess(3))])
disp(['valor teorico ',num2str(n^2*pi^2/(8*a^2))])

C = trapz(x,y.^2);
y = y/(C.^(1/2));
trapz(x,y.^2)
plot(x,y)