function derivadas = e_2017_2018_exer_2_recurso_function(t,solucao,a,b,c)
    derivadas = zeros(2,1);
    X = solucao(1);
    Y = solucao(2);
    derivadas(1) = X*(1-X)-a*X*Y/(X+Y);
    derivadas(2) = b*X*Y/(X+Y)-c*Y;
end