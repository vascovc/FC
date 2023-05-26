clear all

K = 1;m =1;x0 = 1;v0=0;
t0 = 0;
tf = 50;
h = 0.01;
w = sqrt(K/m);

t = t0:h:tf;
N = length(t);

v = zeros(1,N); x = zeros(1,N);
v(1) = v0;
x(1) = x0;

%Euler implicito
A = [ 1 -h/2; w^2*h/2 1];
b = [ x0+h/2*v0 ; v0 - w.^2*h/2*x0 ];
%AZ = b
%determinasse Z com linsolve
for k = 1:(N-1)
    Z = linsolve(A,b);
    x(k+1) = Z(1);
    v(k+1) = Z(2);
    %atualizar b com a nova posicao k+1
    b = [Z(1)+h/2*Z(2); Z(2)-w.^2*h/2*Z(1)];
end

x_exact = x0*cos(w.*t);
v_exact = -x0.*w.*sin(w.*t);
%comparar
plot(t,x_exact,t,x)
plot(t,v_exact,t,v)
%velocidade em func da pos
plot(x,v)

Ec_exact = 1/2 * m *v_exact.^2;
Ep_exact = 1/2 * K *x_exact.^2;
Em_exact = Ec_exact + Ep_exact;

Ec = 1/2 * m *v.^2;
Ep = 1/2 * K *x.^2;
Em = Ec + Ep;
%energia em func do tempo
plot(t,Em_exact,t,Em)