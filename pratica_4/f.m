function derivadas = f(t,solucao,epsilon)
derivadas = zeros(2,1);

y = solucao(1);
v = solucao(2);

derivadas(1) = v;
%derivada de x

derivadas(2) = -epsilon*(y^2-1)*v-y;
%derivade de v
end
