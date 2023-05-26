clear all
close all

h = 0.01;
y0 = 0;
v0 = 2;
F0 = 1.5;
t0=0;tf=100;
epsilon = 0.3;

t = t0:h:tf;
N = length(t);
y = nan(1,N);
v = nan(1,N);
y(1) = y0;v(1) = v0;

fy = @(V) V; %derivada y
fv = @(T,Y,V) F0*cos(1.7*T)-epsilon*(Y^2-1)*V-Y; %derivada velocidade

for k=1:N-1
    r1y = fy( v(k) );
    r1v = fv( t(k), y(k),v(k) );
    
    r2y = fy( v(k) + r1v*h/3 );
    r2v = fv( t(k) + h/3 , y(k) + r1y*h/3, v(k)+r1v*h/3);
    
    r3y = fy( v(k) - r1v*h/3 + r2v*h);
    r3v = fv( t(k) + 2*h/3 , y(k) - r1y*h/3 + r2y*h ,v(k)-r1v*h/3 + r2v*h );
    
    r4y = fy( v(k) + r1v*h -r2v*h + r3v*h );
    r4v = fv( t(k) + h,y(k) + r1y*h-r2y*h + r3y*h,v(k) + r1v*h-r2v*h + r3v*h);
    
    y(k+1) = y(k)+h/8*(r1y + 3*r2y + 3*r3y + r4y);
    v(k+1) = v(k)+h/8*(r1v + 3*r2v + 3*r3v + r4v);
end

subplot(1,2,1)
plot(t,y)
subplot(1,2,2)
plot(y,v)