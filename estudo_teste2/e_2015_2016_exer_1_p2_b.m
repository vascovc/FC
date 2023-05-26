clear all
close all

L = 3; T = 5.0E4; w = 1E5;
x0 = 0;
h = 0.01;
alfa = 5E-8;
x = x0:h:L;
N = length(x);
v = nan(1,N);
y = nan(1,N);
y(1) = 0;
v(1) = -0.01;
B = 0;
guess = [-0.01 0.01];
result = nan(1,2);

for i = 1:1:2
    for k = 1:N-1
        v(k+1) = v(k) + ( 2*alfa*T*guess(i)+alfa*w*x(k)*( L-x(k) ) )*h;
        y(k+1) = y(k) + v(k) * h;
    end
    result(i) = y(N);
end
m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
guess(3) = guess(2) + (B-result(2))/m;
while abs( guess(2)-guess(1) ) > 1E-5
    %para manter o ciclo limpo
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    
    %Euler
    for k = 1:N-1
        v(k+1) = v(k) + ( 2*alfa*T*guess(2)+alfa*w*x(k)*( L-x(k) ) )*h;
        y(k+1) = y(k) + v(k) * h;
    end
    result(2) = y(N);
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    %b = result(2)-m*guess(2); % nao e necessario
    guess(3) = guess(2) + (B-result(2))/m;
    %seria de recomendar usar um valor limite para ver a covergencia
end
subplot(1,2,1)
plot(x,y)
subplot(1,2,2)
plot(x,v)

%Runge-Kutta
clear all
close all

L = 3; T = 5.0E4; w = 1E5;
x0 = 0;
h = 0.01;
alfa = 5E-8;
x = x0:h:L;
N = length(x);
v = nan(1,N);
y = nan(1,N);
y(1) = 0;
v(1) = -0.01;
B = 0;
guess = [-0.01 0.01];
result = nan(1,2);

fy = @(V) V; %derivada posicao
fv = @(y,x) 2*alfa*T*y+alfa*w*x*( L-x ); %derivada tempo
for i = 1:1:2
    for k = 1:N-1
        r1v = fv( guess(i), x(k) );
        r1y = fy( v(k) );
        
        r2v = fv( guess(i),x(k) + r1y*h/2);
        r2y = fy( v(k) + r1v*h/2);
        
        y(k+1) = y(k) + r2y*h;
        v(k+1) = v(k) + r2v*h;
    end
    result(i) = y(N);
end
m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
guess(3) = guess(2) + (B-result(2))/m;
while abs( guess(2)-guess(1) ) > 1E-5
    %para manter o ciclo limpo
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    
    %Euler
    for k = 1:N-1
        r1v = fv( guess(2), x(k) );
        r1y = fy( v(k) );
        
        r2v = fv( guess(2),x(k) + r1y*h/2);
        r2y = fy( v(k) + r1v*h/2);
        
        y(k+1) = y(k) + r2y*h;
        v(k+1) = v(k) + r2v*h;
    end
    result(2) = y(N);
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    %b = result(2)-m*guess(2); % nao e necessario
    guess(3) = guess(2) + (B-result(2))/m;
    %seria de recomendar usar um valor limite para ver a covergencia
end
subplot(1,2,1)
plot(x,y)
subplot(1,2,2)
plot(x,v)
