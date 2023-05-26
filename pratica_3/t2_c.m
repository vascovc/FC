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

fx = @(V) V; %derivada posicao
fv = @(X) -K*X/m; %derivada tempo
for k=1:N-1
    r1v = fv( x(k) );
    r1x = fx( v(k) );
    
    r2v = fv( x(k) + r1x*h/2);
    r2x = fx( v(k) + r1v*h/2);
    
    r3v = fv( x(k) + r2x*h/2);
    r3x = fx( v(k) + r2v*h/2);
    
    r4v = fv( x(k) + r3x*h);
    r4x = fx( v(k) + r3v*h);
    
    v(k+1) = v(k) + ( r1v+2*r2v+2*r3v+r4v )*h/6;
    x(k+1) = x(k) + ( r1x+2*r2x+2*r3x+r4x )*h/6; 
end

subplot(1,2,1)
plot( t,x,t, x0*cos(w.*t) )
subplot(1,2,2)
plot( t,v,t, -w*x0*sin(w.*t) )
figure(2)
plot(x,v,x0*cos(w.*t),-w*x0*sin(w.*t) )

Ec = 1/2 * m *v.^2;
Ep = 1/2 * K *x.^2;
Em = Ec + Ep;
Ec_exact = 1/2 * m *( -w*x0*sin(w.*t) ).^2;
Ep_exact = 1/2 * K *( x0*cos(w.*t) ).^2;
Em_exact = Ec_exact + Ep_exact;
figure(3)
plot(t,Em,t,Em_exact)