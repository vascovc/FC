clear all
close all

m = 0.5; I = 10E-4; K = 5; delta = 10E-3; epsilon = 10E-2;

tf = 100;
h = 0.01;
z0 = 0.1; v0 = 0; teta0 = 0; w0 = 0;
t = 0:h:tf;
N = length(t);
z = nan(1,N);v = nan(1,N);teta = nan(1,N);w = nan(1,N);
z(1) = z0; v(1) = v0; teta(1) = teta0;w(1) = w0;

%A = v w z teta 
A = [1 0 h*K/(2*m) h*epsilon/(4*m);
     0 0 1 -h/2;
     0 1 h*epsilon/(4*I) h*delta/(2*I);
     0 -h/2 0 1];
b = [v(1)-h*K/(2*m)*z(1)-h*epsilon/(4*m)*teta(1);
     z(1)+h/2*v(1);
     w(1)-h*epsilon/(4*I)*z(1)-h*delta/(2*I)*teta(1);
     teta(1)+h/2*w(1)];
for k=1:N-1
    Z = linsolve(A,b);
    v(k+1) = Z(1);
    w(k+1) = Z(2);
    z(k+1) = Z(3);
    teta(k+1) = Z(4);
    
    b = [v(k+1)-h*K/(2*m)*z(k+1)-h*epsilon/(4*m)*teta(k+1);
         z(k+1)+h/2*v(k+1);
         w(k+1)-h*epsilon/(4*I)*z(k+1)-h*delta/(2*I)*teta(k+1);
         teta(k+1)+h/2*w(k+1)];
end