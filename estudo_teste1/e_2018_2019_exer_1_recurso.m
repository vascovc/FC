clear all
close all

R = 0.2; g = 9.8;

w = 5; teta0 = pi/4; v0 = 0;
h = 0.01;
t = 0:h:10;
N = length(t);

teta = nan(1,N); v = nan(1,N);
teta(1) = teta0; v(1) = v0;

f_teta = @(V) V;
fv = @(T) (w^2*cos(T)-g/R)*sin(T);
for k=1:N-1
    r1t = f_teta( v(k) );
    r1v = fv( teta(k) );
    
    r2t = f_teta(v(k) + r1v*h/2 );
    r2v = fv( teta(k) + r1t*h/2);
    
    r3t = f_teta( v(k) + r2v*h/2 );
    r3v = fv( teta(k) + r2t*h/2 );
    
    r4t = f_teta( v(k) + r3v*h );
    r4v = fv( teta(k) + r3t*h );
    
    teta(k+1) = teta(k) + h/6*(r1t + 2*r2t + 2*r3t + r4t);
    v(k+1) = v(k) + h/6*(r1v + 2*r2v + 2*r3v + r4v);
end

plot(t,teta)
for k=2:N-1
    if teta(k-1)<teta(k) && teta(k)>teta(k+1)
        break
    end
end
periodo = interp1( [teta(k-1) teta(k)],[t(k-1) t(k)],teta(k))

