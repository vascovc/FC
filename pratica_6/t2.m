clear all
close all

L = 50; k = 0.93; c = 0.094; ro = 8.9;
T0 = 100;

delta_x = 0.5;
delta_t = 0.1;
tempo_f = 500;
tempo = 0:delta_t:tempo_f;
N_tempo = length(tempo);

eta = k*delta_t/(c*ro*delta_x^2);

x   = 0:delta_x:L;
N_x = length(x);

T = nan(N_x,N_tempo);
% colunas com temperatura de um dado instante
% linhas de uma dada posicao ao longo do tempo

T(:,1) = 100;
T(1,:) = 0;
T(N_x,:) = 0;

% como nao ha extremos sao menos dois pontos
A1 = diag( repmat( 2/eta+2 , 1 , N_x-2 ));
A2 = diag( repmat( -1, 1, N_x-3 ),1);
A3 = diag( repmat( -1, 1, N_x-3 ),-1);
A = A1+A2+A3;

B = nan(N_x-2,1);
for n = 1:N_tempo-1
    for i = 2:N_x-1
        B(i-1) = T(i-1,n)+(2/eta-2)*T(i,n)+T(i+1,n);
    end
    
    %% crank-nicolson
    Z = linsolve(A,B);
    T( 2:(N_x-1),n+1) = Z;
    
    %% sol_sist_trid
    T( 2:(N_x-1),n+1) = sol_sist_trid(A,B);
    
    %% rotina de fatorizacao LU do MAtlab
    [L,U,P] = lu(A);%mandar para fora do ciclo para eficiencia
    y = L\B;
    T( 2:(N_x-1),n+1) = U\y;
    
end

mesh(tempo,x,T); xlabel('tempo'); ylabel('x');zlabel('temperatura');
figure(2)
contourf(x,tempo,T')
colorbar