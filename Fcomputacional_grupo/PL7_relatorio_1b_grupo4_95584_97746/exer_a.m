%PL7
%Tiago Santos 95584
%Vasco Costa  97746
clear all
close all

L = 4; 
h = 0.025; % espacamento
N_max_iter = 1E6;% valor maximo de iteracoes, evitar calculos infinitos
tol = 1E-7;%tolerancia

y = -L/2:h:L/2;
x = -L/2:h:L/2;

N = length(x);%tamanho segundo uma direcao vai ser o M
z_old = zeros(N,N);

z_old(1,:) = 10;
z_old(end,:) = 12;
z_old(:,1) = 11+2/L.*x;
z_old(:,end) = 11+2/L.*x;

f = zeros(N,N);%faz se o calculo de todos os f, dentro da circunferencia vai ser 1 e fora vai valer 0, evita if
for j = 1:N
    for i = 1:N
        if (y(i)^2+ x(j)^2) < (L^2/9)
            f(j,i) = 1;
        end
    end
end

z_inicial = z_old;

figure(1)
meshc(x,y,z_old),xlabel('x'),ylabel('y'),title('condicoes iniciais')

%% Jacobi
z_new = z_old;
for iter = 1:N_max_iter
    for i = 2:(N-1)
        for j = 2:(N-1)
            z_new(j,i) = 1/4*( z_old(j+1,i) + z_old(j-1,i) + z_old(j,i+1) + z_old(j,i-1) - h^2*f(j,i) );
        end
    end
    num = sqrt(sum(sum((z_new-z_old).^2)));
    den = sqrt(sum(sum(z_new.^2)));
    if (num/den) < tol
        break
    end
    z_old = z_new;
end
disp(['Jacobi numero de iter:',num2str(iter)])
figure(2)
meshc(x,y,z_new),xlabel('x'),ylabel('y'),title('Método de Jacobi')


z_old = z_inicial;
%% Gauss-Seidel
z_new = z_old;
for iter = 1:N_max_iter
    for i = 2:(N-1)
        for j = 2:(N-1)              
            z_new(j,i) = 1/4*( z_new(j+1,i) + z_new(j-1,i) + z_new(j,i+1) + z_new(j,i-1) - h^2*f(j,i) );       
        end
    end    
    num = sqrt(sum(sum((z_new-z_old).^2)));
    den = sqrt(sum(sum(z_new.^2)));
    if (num/den) < tol
        break
    end
    z_old = z_new;
end
disp(['Gauss-Seidel numero de iter:',num2str(iter)])

figure(3)
meshc(x,y,z_new),xlabel('x'),ylabel('y'),title('Método de Gauss-Seidel')


z_old = z_inicial;
%% sobre-relaxacao sucessiva
alfa = 2-2*pi/N;
z_new = z_old;
for iter = 1:N_max_iter
    for i = 2:(N-1)
        for j = 2:(N-1)
            z_new(j,i) = (1-alfa)*z_old(j,i)+alfa/4*( z_new(j+1,i) + z_new(j-1,i) + z_new(j,i+1) + z_new(j,i-1) - h^2*f(j,i) );      
        end
    end
    num = sqrt(sum(sum((z_new-z_old).^2)));
    den = sqrt(sum(sum(z_new.^2)));
    if (num/den) < tol
        break
    end
    z_old = z_new;
end
disp(['sobre-relaxacao iter:',num2str(iter)])
figure(4)
meshc(x,y,z_new),xlabel('x'),ylabel('y'),title('Método da sobre-relaxação sucessiva')
