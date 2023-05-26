function derivadas = f(t,solucao,m,K,alfa)
derivadas = zeros(2,1);

x = solucao(1);
v = solucao(2);

derivadas(1) = v;
%derivada de x

derivadas(2) = -K/m*x*(1+3/2*alfa*x);
%derivade de v
end