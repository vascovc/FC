clear all
close all

h = 0.1;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);
t = t0:h:tf;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = x0;v(1) = v0;

for k=1:N-1
    r1x = v(k);
    r1v = -w^2 * x(k);
    
    r2x = v(k) + r1v*h/2;
    r2v = -w^2 * ( x(k) + r1x*h/2 );
    
    v(k+1) = v(k) + r2v*h;
    x(k+1) = x(k) + r2x*h;
end

subplot(1,2,1)
plot( t,x,t, x0*cos(w.*t) )
subplot(1,2,2)
plot( t,v,t, -w*x0*sin(w.*t) )

