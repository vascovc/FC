clear all

K = 1;m =1;x0 = 1;v0 = 0;
t0 = 0;
tf = 50;
h = 0.001;

t = t0:h:tf;
N = length(t);

v = zeros(1,N); x = zeros(1,N);
v(1) = v0;
x(1) = x0;

for k = 1:(N-1)%Euler-Cromer
    v(k+1) = v(k)-(K/m)*x(k)*h;
    x(k+1) = x(k)+v(k+1)*h;%%Euler-cromer so faz com que aqui se meta mais 1
end

w = sqrt(K/m);

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