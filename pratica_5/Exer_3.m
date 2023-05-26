clear all
close all

R = 1E-3; Q = 2.1E6; lambda = 0.1;
h = 1E-5;
r = 0:h:R;
N = length(r);
A1 = diag( repmat(-2,[1 N]) );
A1(end,end) = 1;
A2 = diag( 1+h./(2.*r(1:end-1)) ,1);
A2(1,2) = 2;
A3 = diag( 1-h./(2.*r(2:end)) ,-1);
A3(end,end-1) = 0;
A = A1+A2+A3;

b = repmat(-h^2*Q/lambda,[N 1]);
b(end) = 20;

T = linsolve(A,b);
plot(r,T)
