function derivadas = e_2015_2016_exer_1_func(t,solucao,alfa,T,w,x,L)
    derivadas = zeros(2,1);
    y = solucao(1);
    v = solucao(2);
    derivadas(1) = v;
    %derivada de x e a velocidade

    derivadas(2) = 2*alfa*T*y+alfa*w*x*(L-x);
end