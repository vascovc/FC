clear all
close all

L = 1;
h = 0.025;
N_max_iter = 1E4;
tol = 1E-7;
delta_x = h;
delta_y = h;

y = -L:delta_y:L;
x = -L:delta_x:L;

%[X,Y] = meshgrid(x,y);
N = length(x);
alfa = 1.65:0.01:1.96;
alfa(end+1) = 2/(1+pi/N);
iteracoes = nan(1,length(alfa));

for ind_alfa = 1:length(alfa)
    V_old = zeros(N,N);
    ind_min = L/(2*h)+1;
    ind_max = 3*L/(2*h)+1;

    V_old( ind_min:ind_max , ind_min:ind_max ) = 1;

    cond_fronteira = V_old == 0;

    V_new = V_old;

    for iter = 1:N_max_iter
        for i = 2:(length(x)-1)
            for j = 2:(length(y)-1)
                if true(cond_fronteira(j,i))     
                    V_new(j,i) = (1-alfa(ind_alfa))*V_old(j,i)+alfa(ind_alfa)/4*( V_new(j+1,i) + V_new(j-1,i) + V_new(j,i+1) + V_new(j,i-1) );
                end        
            end
        end
        num = sqrt(sum(sum((V_new-V_old).^2)));
        den = sqrt(sum(sum(V_new.^2)));
        if (num/den) < tol
            break
        end
        V_old = V_new;
    end
    iteracoes(ind_alfa) = iter;
end
plot(alfa(1:end-1),iteracoes(1:end-1),'.',alfa(end),iteracoes(end),'*')

[u,v] = min(iteracoes);
disp(['menor num de iter ',num2str(iteracoes(v)),' com alfa = ',num2str(alfa(v)) ]);

