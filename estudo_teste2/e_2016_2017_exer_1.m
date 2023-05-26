clear all
close all

L = 5;
gama = -1.5;
epsilon = 2;
h = 0.01;
x = 0:h:L;
N=length(x);
A1 = diag( repmat(-4*epsilon,[1 N]) );
A1(1,1) = 1;
A1(end,end) = 1;
A2 = diag(repmat(2*epsilon+h,[1 N-1]),-1);
A2(end,end-1) = 0; 
A3 = diag(repmat(2*epsilon-h,[1 N-1]),1);
A3(1,2) = 0;
A = A1+A2+A3;

B = nan(N,1);
B(:,1) = 2*h^2*gama;
B(1) = 20;
B(end) = 50;

T = linsolve(A,B);
plot(x,T)
hold on
plot(x,B(1)-gama*x+(B(end)-B(1)+gama*L)*(exp(x/epsilon)-1)/(exp(L/epsilon)-1))

