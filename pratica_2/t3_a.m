clear all

alfa = -0.1; m=1; 
v0 = 1;x0 = 1;
K = 1;
t0 = 0;
tf = 50;
h = 0.01;

t = t0:h:tf;
N = length(t);
v = nan(N,1); x = nan(N,1);
v(1) = v0;
x(1) = x0;
%Euler Cromer
for k = 1:(N-1)
    a = -K/m *(x(k) +2*alfa*x(k)^3);
    v(k+1) = v(k) + a*h;
    x(k+1) = x(k) + v(k+1)*h;
end
contador = 0;
for k=2:N-1
    if x(k-1) <= x(k) && x(k) >= x(k+1)
        contador = contador+1;
        aux = lagr( t(k-1:k+1), x(k-1:k+1) );
        t_max(contador) = aux(1);
        x_max(contador) = aux(2);
    end
end

A = mean(x_max)
P = mean(diff(t_max))