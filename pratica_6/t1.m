clear all
close all

L = 50; k = 0.93; c = 0.094; ro = 8.9;
T0 = 100;

delta_x = 0.5;
delta_t = 0.1;
tempo_f = 500;
tempo = 0:delta_t:tempo_f;
N_tempo = length(tempo);

x   = 0:delta_x:L;
N_x = length(x);

T = nan(N_x,N_tempo);
% colunas com temperatura de um dado instante
% linhas de uma dada posicao ao longo do tempo

T(:,1) = 100;
T(1,:) = 0;
T(N_x,:) = 0;
%[ 0  0  0 ...
% 100 
% 100 
%  .
%  .
% 100
%  0  0  0 ...]
%   T(posicao_x,tempo)
for n = 1:N_tempo-1 % o tempo comeca logo na primeira coluna
    for i = 2:N_x-1 % a temperatura calcula se so a partir do 2
        T(i,n+1) = T(i,n) + k*delta_t/(c*ro*delta_x^2) * ( T(i-1,n) - 2*T(i,n) + T(i+1,n) );
    end 
end

mesh(tempo,x,T); xlabel('tempo'); ylabel('x');zlabel('temperatura');
figure(2)
contourf(x,tempo,T')
colorbar

%%c
%historia de encontrar os pontos que estao mais perto de L/4
index_x = find( abs( ( x - L/4 ) ) < delta_x/2 );

figure(3)
plot(tempo,T(index_x,:)),xlabel('tempo'),ylabel('temperatura')

[~,index_tempo] = min( abs( T(index_x,:)-50 ) );
%o ponto mais perto da temperatura menos 50 igual a 0, ou seja o mais perto
%de 50
disp(['tempo de 50º ',num2str(tempo(index_tempo)),' segundos']);