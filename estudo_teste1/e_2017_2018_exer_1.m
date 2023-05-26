clear all
close all

w = 4;
teta0 = 1.1;T = 200;

a=T*sin(teta0)/w;

h=0.01;
x = 0:h:100;
N = length(x);
x = nan(1,N);
y = nan(1,N);
teta = nan(1,N);
teta(1) = teta0;x(1) = 0;y(1)=0;

fy = @(TETA) -cot(TETA);
fv = @(V) 1/a*sqrt(1+V^2);
for k=1:N-1
    r1y = fy(teta(k));
    r1v = fv(y(k));
    
    r2y = fy(teta(k)+r1y*h/2);
    r2v = fv(y(k)+r1v*h/2);
    
    r3y = fy(teta(k)+r2y*h/2);
    r3v = fv(y(k)+r2v*h/2);
    
    r4y = fy(teta(k)+r3y*h);
    r4v = fv(y(k)+r3v*h);
    
    y(k+1) = y(k) +h/6*(r1y+2*r2y+2*r3y+r4y);
end