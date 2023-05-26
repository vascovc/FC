clear all
close all

N = 8;
i = 1:1:N;

%metodo que deu mais trabalho, dava definir diagonais
a = [-1 0 0 sqrt(2)/2 1 0 0 0;
    0 -1 0 sqrt(2)/2 0 0 0 0;
    0 0 -1 0 0 0 1/2 0;
    0 0 0 -sqrt(2)/2 0 -1 1/2 0;
    0 0 0 0 -1 0 0 1;
    0 0 0 0 0 1 0 0;
    0 0 0 -sqrt(2)/2 0 0 sqrt(3)/2 0;
    0 0 0 0 0 0 -sqrt(3)/2 -1];
b = zeros(N,1);
b(6) = 1E4;

x_old = ones(N,1);
x_new = x_old;

kmax = 500;
for k=1:kmax
    for i = 1:N
        soma = 0;
        for j = 1:N
            if j ~= i
                soma = soma + a(i,j)*x_old(j);
            end
        end
        x_new(i) = -1/a(i,i) * soma + 1/a(i,i) * b(i);
    %tambem dava para fazer em dois ciclos for de 1 ate i-1 e de 1+1 ate N
    end
    
    if max(abs(x_new-x_old))/max(abs(x_new)) < 1E-7
        break
    end
    x_old = x_new;
end
disp(['numero de iteracoes ',num2str(k)])