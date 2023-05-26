clear all
close all

L = 50; k = 0.93; ro = 8.9;
T0 = 100;

delta_x = 0.5;
delta_t = 0.1;
tempo_f = 250;
tempo = 0:delta_t:tempo_f;
N_tempo = length(tempo);

x   = 0:delta_x:L;
N_x = length(x);

c = nan(N_x,1);
c(:,1) = 0.094;
c( floor(N_x/2):end ) = 0.188;

eta = k.*delta_t./(c.*ro.*delta_x^2);
T = nan(N_x,N_tempo);
% colunas com temperatura de um dado instante
% linhas de uma dada posicao ao longo do tempo

T(:,1) = 50*sin(2*pi*x/L);
T(1,:) = 0;
T(N_x,:) = 0;

A1 = diag( 2./eta(2:end-1) + 2);
A2 = diag( repmat( -1, 1, N_x-3 ),1);
A3 = diag( repmat( -1, 1, N_x-3 ),-1);
A = A1+A2+A3;

B = nan(N_x-2,1);
for n = 1:N_tempo-1
    for i = 2:N_x-1
        B(i-1) = T(i-1,n)+( 2/eta(i) - 2)*T(i,n) + T(i+1,n);
    end
    Z = linsolve(A,B);
    T( 2:(N_x-1),n+1) = Z;
end

mesh(tempo,x,T); xlabel('tempo'); ylabel('x');zlabel('temperatura');
figure(2)
contourf(x,tempo,T')
colorbar