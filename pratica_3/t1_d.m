clear all
close all

h = 0.01;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);
t = t0:h:tf;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = x0;v(1) = v0;

fx = @(T,X,V) V; %derivada posicao
fv = @(T,X,V) -K*X/m; %derivada tempo
for k=1:N-1
    r1x = fx( t(k),x(k),v(k) );
    r1v = fv( t(k),x(k),v(k) );
    
    r2x = fx( t(k) +h/2 ,x(k) + r1x * h/2 ,v(k) + r1v*h/2 );
    r2v = fv( t(k) +h/2 ,x(k) + r1x * h/2 ,v(k) + r1v*h/2 );
    
    v(k+1) = v(k) + r2v*h;
    x(k+1) = x(k) + r2x*h;
end

subplot(1,2,1)
plot( t,x,t, x0*cos(w.*t) )
subplot(1,2,2)
plot( t,v,t, -w*x0*sin(w.*t) )

%e

clear all
close all

h = 0.01;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);
t = t0:h:tf;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = x0;v(1) = v0;

fx = @(V) V; %derivada posicao
fv = @(X) -K*X/m; %derivada tempo
for k=1:N-1
    r1x = fx( v(k) );
    r1v = fv( x(k) );
    
    r2x = fx( v(k) + r1v * h/2 );
    r2v = fv( x(k) + r1x * h/2 );
    
    v(k+1) = v(k) + r2v * h;
    x(k+1) = x(k) + r2x * h;
end

subplot(2,2,1)
plot( t,x,t, x0*cos(w.*t) )
subplot(2,2,2)
plot( t,v,t, -w*x0*sin(w.*t) )

subplot(2,2,3)
plot(x,v,x0*cos(w.*t),-w*x0*sin(w.*t))