function derivadas = f(t,solucao)
derivadas = zeros(2,1);
% Esta linha � necess�ria e faz do vetor de sa�da um vetor coluna.
K=16;
m=1;

% O vetor solucao tem os valores de x e v para o tempo
% t em que a fun��o � chamada pela rotina ode45.
derivadas(1) = solucao(2);
%derivada de x e a velocidade

% Escreva acima a express�o da derivada de x em fun��o
%de solucao(1) e de solucao(2).
derivadas(2) = -K/m * solucao(1);
%a derivada da velocidade e -k sobre m vezes x
% Escreva a acima express�o da derivada de v em fun��o de solucao(1) e de solucao(2).
end
