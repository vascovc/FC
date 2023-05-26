clear all
close all

tf = 0.5; z0 = 50E-2;A = 8.3;B = 21;
h = 0.05;
t = 0:h:tf;
N = length(t);

z = nan(1,N);z(1) = z0;
v = nan(1,N);v(1) = 0;

m_a = [1 -h/2; 0 B*h/2+1];
m_b = [z(1)+v(1)*h/2;-A*h+v(1)*(1-B*h/2)];
for k=1:N-1
    Z = linsolve(m_a,m_b);
    z(k+1) = Z(1);
    v(k+1) = Z(2);
    m_b = [z(k+1)+v(k+1)*h/2;-A*h+v(k+1)*(1-B*h/2)];
end

subplot(1,2,1)
plot(t,v,t,A/B*(exp(1).^(-B.*t)-1))