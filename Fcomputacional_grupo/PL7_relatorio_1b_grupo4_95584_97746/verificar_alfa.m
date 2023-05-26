%PL7
%Tiago Santos 95584
%Vasco Costa  97746
clear all
close all

L = 4;
h = 0.025;
N_max_iter = 1E6;
tol = 1E-7;

y = -L/2:h:L/2;
x = -L/2:h:L/2;

N = length(x);
z_old = zeros(N,N);

z_old(1,:) = 10;
z_old(end,:) = 12;
z_old(:,1) = 11+2/L.*x;
z_old(:,end) = 11+2/L.*x;

f = zeros(N,N);
for j = 1:N
    for i = 1:N
        if (y(i)^2+ x(j)^2) < (L^2/9)
            f(j,i) = 1;
        end
    end
end
z_inicial = z_old;

alfa = 1.80:0.001:1.99;
num_total_it = nan(1,length(alfa));
for ind_alfa = 1:length(alfa)
    z_old = z_inicial;
    z_new = z_old;
    for iter = 1:N_max_iter
        for i = 2:(N-1)
            for j = 2:(N-1)
                z_new(j,i) = (1-alfa(ind_alfa))*z_old(j,i)+alfa(ind_alfa)/4*( z_new(j+1,i) + z_new(j-1,i) + z_new(j,i+1) + z_new(j,i-1) - h^2*f(j,i) );      
            end
        end
        num = sqrt(sum(sum((z_new-z_old).^2)));
        den = sqrt(sum(sum(z_new.^2)));
        if (num/den) < tol
            break
        end
        z_old = z_new;
    end
    num_total_it(ind_alfa) = iter;
    disp(['alfa = ',num2str(alfa(ind_alfa))])
end
plot(alfa,num_total_it,'.')

% para ver o alfa_opt
%alfa = 2/(1+pi/N); da aula pratica
alfa = 2-2*pi/N;
z_old = z_inicial;
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

hold on
plot(alfa,iter,'*')
xlabel('alfa'),ylabel('iteracoes'),title('iteracoes em funcao do valor de alfa')