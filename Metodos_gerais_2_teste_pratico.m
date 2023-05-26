%% Metodo do shooting Folha 5 - 24
m = ( result(2)-result(1) )/( guess(2) -guess(1) );
b = result(2) - m * guess(2);
% B vai ser o valor que nos queremos aproximar
guess(3) = guess(2) + ( B - result(2) )/m;
%fazer disto ciclo

%Folha 5 - Exer_1_3.m
%Folha de revisoes - exer1.m
%% Diferencas finitas Linsolve Folha 5 - 8
y' = ( y(x+h) - y(x) )/h +O(h); %diferenca avancada/progressiva
y' = ( y(x) - y(x-h) )/h +O(h); %diferenca retardada/regressiva
y' = ( y(x+h) - y(x-h) )/(2*h) +O(h^2); %diferenca centrada

y'' = ( y(x+h) -2*y(x) + y(x-h) )/(h^2)+O(h^2); %diferenca centrada

%substituir, multiplicar por h^2 se tiver a segunda derivada
%organizar a matriz de forma a se ter os n-1 n n+1 e depois ver

%condicoes fronteira
%Sabendo a derivada
T'(1) = 0 % y' = ( y(x+h) -y(x) )/h
0 = (T(2) - T(1))/h => -T1 + T2 = 0
%sabendo o valor
T(N) = 20
%construir a matriz A[1 N] e depois a b, obter
%T com linsolve

%Folha 5 - exercicio 3.m
%% Valores proprios Folha 5 - 28
%ver que da para encontrar algo comum a todo os y(k)
%obter a matriz A [2 N-1] tem ordem [x-h x x+h]
sol = eigs(A,3,'sm'); %3 valores proprios

%Folha 5 - exercicio 2.m
%% Calor Euler
%Folha 6 - t1.m
%% Crank-Nicolson
  %linsolve
    Z = linsolve(A,B);
    T( 2:(N_x-1),n+1) = Z;
  %sol_sist_trid
    T( 2:(N_x-1),n+1) = sol_sist_trid(A,B);
  %fatorizacao lu
    [L,U,P] = lu(A);%mandar para fora do ciclo para eficiencia
    y = L\B;
    T( 2:(N_x-1),n+1) = U\y;
    
%Folha 6 - t2.m
%Folha 6 - t3_.m
%Folha revisoes - exer3.m %aproximacao paraxial
%% Transformadas de Fourier
delta_w = 2*pi/(N*delta_t);
t = 0:delta_t:(N-1)*delta_t;
delta_t = 2*pi/(N*delta_w);
%densidade espetral - F(k)^2
dens_espt = (delta_t*abs(Z)).^2;

%freq_nyquist
w_max = N*delta_w/2 = pi/delta_t;
f_max = 1/(2*delta_t);

%aliasing, para n ocorrer, aumentar w_max => diminuir delta_t
delta_t < 1/(2*f_high) = T_low/2;

% sem shift
w_max = (N-1)*delta_w;
w_min = 0;
w(:,1) = w_min:delta_w:w_max;
% com shift para n par
w_max = (N/2-1)*delta_w;
w_min = (-N/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;
%com shift para n impar adiciona-se delta_w a ambos
w_max = (N/2-1/2)*delta_w;
w_min = (-N/2+1/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;

Z = fftshift( fft(y) );

%Folha 7 exer1_a.m / exer1_b.m / exer1_c.m
%Folha 7 exer2.m (com uso da inversa para retomar ao normal)
%Folha de revisoes exer2.m