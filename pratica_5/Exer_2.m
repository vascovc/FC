clear all
close all

u = 10^-3; L = 1; T = 10^3;
n = [1 2 3];
wn = n*pi/L*sqrt(T/u);
x0 = 0;
h = 0.01;
x = x0:h:L;
N = length(x);

A1 = diag( repmat(-2,1,N-2) );
A2 = diag( ones(1,N-3),1 );%linha superior ao menos 2
A3 = diag( ones(1,N-3),-1 );%linha inferior ao menos 2

A = A1 + A2 + A3;%matriz A
sol = eigs(A,3,'sm');%3 para obter os primeiros 3 modos de vibracao
sol = sqrt( -sol*T/(u*h^2) );

[vec,val] = eigs(A,3,'sm');%para se terem os pontos do y em vec
figure(1)
plot(x(2:end-1),vec(:,1:3))%imprime os 3 modos