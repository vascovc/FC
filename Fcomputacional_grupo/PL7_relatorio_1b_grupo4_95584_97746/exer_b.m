%PL7
%Tiago Santos 95584
%Vasco Costa  97746
% contar o tempo
% tic - comeca
% tempo = toc - conta o fim do tempo

clear all
close all

L = 4;
N_max_iter = 1E6;
tol = 1E-6;

h = zeros(1,30);
for ind = 1:30
    h(ind) = (L/2)/( ind*5+60 );%(L/2)/80 da os 0.025 e assim temos uma serie de h que sao multiplos de h/2
end
%h = [0.250 0.100 0.050 0.025 0.010 0.0050 0.0025];

N_size = nan(1,length(h));%para os tamanho de M
time = nan(1,length(h));%guardar os tempos
iteracoes = nan(1,length(h));%guardar as iteracoes
for ind_h = 1:length(h)
    tic%comecamos aqui o tempo

    y = -L/2:h(ind_h):L/2;
    x = -L/2:h(ind_h):L/2;

    N = length(x);
    N_size(ind_h) = N;
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

    z_new = z_old;
    for iter = 1:N_max_iter
        for i = 2:(N-1)
            for j = 2:(N-1)               
                z_new(j,i) = 1/4*( z_new(j+1,i) + z_new(j-1,i) + z_new(j,i+1) + z_new(j,i-1) - h(ind_h)^2*f(j,i) );         
            end
        end

        num = sqrt(sum(sum((z_new-z_old).^2)));
        den = sqrt(sum(sum(z_new.^2)));
        if (num/den) < tol
            break
        end
        z_old = z_new;
    end
    time(ind_h) = toc;%acaba aqui o tempo
    iteracoes(ind_h) = iter;
    disp(['iteracao: ',num2str(ind_h)])%para ajudar a ver onde estamos
end
figure(1)
plot( log(N_size),log(iteracoes),'o' ),xlabel('ln(M)'),ylabel('ln(iteracoes)'),title('iteracoes')
lsline%reta de melhor ajuste
p = polyfit( log(N_size),log(iteracoes),1 );
display(['declive iteracoes = ',num2str(p(1))])

figure(2)
plot( log(N_size),log(time),'o' ),xlabel('ln(M)'),ylabel('ln(time)'),title('tempo')
lsline
p = polyfit( log(N_size),log(time),1 );
display(['declive tempo = ',num2str(p(1))])