clear all
close all

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
v(1) = 0.2;

for k = 1:N-1
    v(k+1) = v(k) + -wn(1)^2*u/T*y(k) * h;
    y(k+1) = y(k) + v(k+1) * h;
end
subplot(1,2,1)
plot(x,y)
%b
wn(1) = 3500;
for k = 1:N-1
    v(k+1) = v(k) + -wn(1)^2*u/T*y(k) * h;
    y(k+1) = y(k) + v(k+1) * h;
end
subplot(1,2,2)
plot(x,y)
