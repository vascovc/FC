clear all
close all

h = 0.01;
y0 = 2;
v0 = 7;
t0=0;tf=100;
epsilon = 1;

t = t0:h:tf;
N = length(t);
y = nan(1,N);
v = nan(1,N);
y(1) = y0;v(1) = v0;

fy = @(V) V; %derivada y
fv = @(Y,V) -epsilon*(Y^2-1)*V-Y; %derivada velocidade

for k=1:N-1
    r1y = fy( v(k) );
    r1v = fv( y(k),v(k) );
    
    r2y = fy( v(k) + r1v * h);
    r2v = fv( y(k) + r1y * h ,v(k) + r1v * h);
    
    r3y = fy( v(k) + r1v *h/4 + r2v *h/4);
    r3v = fv( y(k) + r1y *h/4 + r2y *h/4,v(k) + r1v *h/4 + r2v *h/4);
    
    y(k+1) = y(k) + h/6*( r1y + r2y + 4*r3y);
    v(k+1) = v(k) + h/6*( r1v + r2v + 4*r3v);
end

subplot(1,2,1)
plot(t,y)
subplot(1,2,2)
plot(y,v)