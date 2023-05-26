clear all
close all
clc

K = 2.5; m = 0.5; x0 = 0.10 ;v0 = 0; h = 0.01;

t = 0:h:10;
N = length(t); v = zeros(1,N); x = zeros(1,N);
v(1) = v0;x(1) = x0;

for k = 1:N-1
    v(k+1) = v(k) + (-K/m)*x(k)*h; %aceleracao altera em cada ciclo
    x(k+1) = x(k) + v(k)*h;
end
w = sqrt(K/m);%frequencia
%x(t) = x0*cos(w*t)
plot(t,x, t, x0*cos(w.*t) )
figure(2)
%v(t) = -x0*w*sin(w*t)
plot(t,v, t, -x0*w*sin(w.*t) )

%o erro vem do overshoot, amplitude aumenta mas periodo mantem
