clear all

K = 1;m =1;x0 =1;v0=0;
t0 = 0;
tf = 50;
h = 0.01;

t = t0:h:tf;
N = length(t);

v = zeros(1,N); x = zeros(1,N);
v(1) = v0;
x(1) = x0;
%Euler-Cromer
for k = 1:(N-1)
    v(k+1) = v(k)-(K/m)*x(k)*h;
    x(k+1) = x(k)+v(k+1)*h;%Euler-cromer
end

w = sqrt(K/m);

contador = 0;
for k=2:(N-1)
    if x(k-1) <= x(k) && x(k) >= x(k+1)
        contador = contador+1;
        aux = lagr( t(k-1:k+1), x(k-1:k+1) );
        t_max(contador) = aux(1);
        x_max(contador) = aux(2);
    end
end

A = mean( x_max )
p = mean( diff(t_max) )%periodo