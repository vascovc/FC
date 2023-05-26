clear all
close all

h = 0.01;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=1000;

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

%euler
x_Euler = nan(1,N);
v_Euler = nan(1,N);
x_Euler(1) = x0;
v_Euler(1) = v0;
for k = 1:(N-1)%Euler
    v_Euler(k+1) = v_Euler(k)-(K/m)*x_Euler(k)*h;
    x_Euler(k+1) = x_Euler(k)+v_Euler(k)*h;
end

subplot(2,2,1)
plot( t,x,t, x0*cos(w.*t),t,x_Euler)
title('posicao')
subplot(2,2,2)
plot( t,v,t, -w*x0*sin(w.*t) ,t,v_Euler)
title('velocidade')

subplot(2,2,3)
plot(x,v,x_Euler,v_Euler,x0*cos(w.*t),-w*x0*sin(w.*t))
title('fase')

Ec = 1/2 * m *v.^2;
Ep = 1/2 * K *x.^2;
Em = Ec + Ep;

Ec_Euler = 1/2 * m *v_Euler.^2;
Ep_Euler = 1/2 * K *x_Euler.^2;
Em_Euler = Ec_Euler + Ep_Euler;

Ec_exact = 1/2 * m *( -w*x0*sin(w.*t) ).^2;
Ep_exact = 1/2 * K *( x0*cos(w.*t) ).^2;
Em_exact = Ec_exact + Ep_exact;

subplot(2,2,4)
plot(t,Em,t,Em_Euler,t,Em_exact)
title('energia')