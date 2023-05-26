clear all
close all

epsilon = 5.0;L =0.1;R =10;C=1E-3;
v0=0;

h=0.00001;
t = 0:h:3*R*C;
N = length(t);

a = 1/(R*C)+R/L;
b = 2/(L*C);
c = epsilon/(L*C);

v = nan(1,N);v(1) = v0;
x = nan(1,N);x(1) = 0;
for k=1:N-1
   v(k+1) = v(k) +(-a*v(k)-b*x(k)+c)*h;
   x(k+1) = x(k) + v(k)*h;
   if abs(x(k)-x(k+1)-epsilon/2) < 10E-6
       break
   end
end
plot(t,x)

%b
epsilon = 5.0;L =0.1;R =10;C=1E-3;
v0=0;

h=0.00001;
t = 0:h:3*R*C;
N = length(t);

a = 1/(R*C)+R/L;
b = 2/(L*C);
c = epsilon/(L*C);

v = nan(1,N);v(1) = v0;
x = nan(1,N);x(1) = 0;

A = [h*b a*h+1;1 -1*h];
B = [v(1)+h*c;x(1)];
for k=1:N-1
   Z = linsolve(A,B);
   v(k+1) = Z(2);
   x(k+1) = Z(1);
   B = [v(k+1)+h*c;x(k+1)];
   if abs(x(k)-x(k+1)-epsilon/2) < 10E-6
       break
   end
end
hold on
plot(t,x)

%c
epsilon = 5.0;L =0.1;R =10;C=1E-3;
v0=0;

h=0.00001;
t = 0:h:3*R*C;
N = length(t);

a = 1/(R*C)+R/L;
b = 2/(L*C);
c = epsilon/(L*C);

v = nan(1,N);v(1) = v0;
x = nan(1,N);x(1) = 0;

fv = @(X,V) -a*V-b*X+c;
fx = @(V) V;
for k=1:N-1
    r1x = fx(v(k));
    r1v = fv(x(k),v(k));
    
    r2x = fx(v(k)+h/2*r1v);
    r2v = fv(x(k)+h/2*r1x,v(k)+h/2*r1v);
    
    x(k+1) = x(k)+r2x*h;
    v(k+1) = v(k)+r2v*h;
end
plot(t,x)