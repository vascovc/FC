clear all
close all
%nao esta certo
L = 2;D = 1;
delta_x = 0.1;
delta_t = 0.01;
tempo_f = 20;
tempo = 0:delta_t:tempo_f;
N_tempo = length(tempo);

x   = -L/2:delta_x:L/2;
N_x = length(x);

ro = nan(N_x,N_tempo);
% colunas com ro de um dado instante
% linhas de uma dada posicao ao longo do tempo
ro(:,1) = 2*cos(pi*x/L);
ro(1,:) = 0;
ro(N_x,:) = 0;

for n = 1:N_tempo-1
    for i = 2:N_x-1
        ro(i,n+1) = ro(i,n) + delta_t/(delta_x^2)*( ro(i+1,n)-2*ro(i,n)+ro(i+1,n) );
    end
end

mesh(tempo,x,ro)
figure(2)
contourf(x,tempo,ro')
colorbar