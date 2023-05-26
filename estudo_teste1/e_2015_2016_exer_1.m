clear all
close all

ro = 1.225; C = 0.508; m = 57E-3;A = 3.5E-4;g = 9.8;
z0 = 0.5; v0 = 15;

h=0.01;
t = 0:h:5;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = z0; v(1) = v0;

fx = @(V) V;
fv = @(V) -g-0.5/m*C*ro*A*norm(V)^2*V;
for k=1:N-1
    r1x = fx( v(k) );
    r1v = fv( v(k) );
    
    r2x = fx( v(k) + r1v*h/2 );
    r2v = fv( v(k) + r1v*h/2 );
    
    r3x = fx( v(k) + r2v*h/2 );
    r3v = fv( v(k) + r2v*h/2 );
    
    r4x = fx( v(k) + r3v*h );
    r4v = fv( v(k) + r3v*h );
    
    x(k+1) = x(k) + h/6*(r1x + 2*r2x + 2*r3x + r4x);
    v(k+1) = v(k) + h/6*(r1v + 2*r2v + 2*r3v + r4v);
end

subplot(1,2,1)
plot(t,x)
subplot(1,2,2)
plot(t,v)

%%b
for k = 2:N-1
    if x(k-1) < x(k) && x(k) > x(k+1)
        break
    end
end
disp(['maior z ',num2str( x(k) )])

%%c
for k=2:N-1
    if(x(k)<0)
        break
    end
end

tempo_impact = interp1( [x(k-1) x(k)], [t(k-1) t(k)],0);
disp(['tempo impacto ',num2str( tempo_impact )] )

%%d
ro_0 = 1.225; C = 0.508; m = 57E-3;A = 3.5E-4;g = 9.8;
z0 = 7000; v0 = 200;

h=0.01;
t = 0:h:2000;
N = length(t);
x = nan(1,N);
v = nan(1,N);
x(1) = z0; v(1) = v0;

ro = @(z) ro_0*exp(-z/z0);
fx = @(V) V;
fv = @(z,V) -g-0.5/m*C*ro(z)*A*norm(V)^2*V;
for k=1:N-1
    r1x = fx( v(k) );
    r1v = fv( x(k),v(k) );
    
    r2x = fx( v(k) + r1v*h/2 );
    r2v = fv( x(k) + r1x*h/2, v(k) + r1v*h/2 );
    
    r3x = fx( v(k) + r2v*h/2 );
    r3v = fv( x(k) + r2x*h/2,v(k) + r2v*h/2 );
    
    r4x = fx( v(k) + r3v*h );
    r4v = fv( x(k) + r1x*h, v(k) + r3v*h );
    
    x(k+1) = x(k) + h/6*(r1x + 2*r2x + 2*r3x + r4x);
    v(k+1) = v(k) + h/6*(r1v + 2*r2v + 2*r3v + r4v);
end
disp( ['max z novo ',num2str(max(x))])
for k=2:N-1
    if(x(k)<0)
        break
    end
end

tempo_impact = interp1( [x(k-1) x(k)], [t(k-1) t(k)],0);
disp(['tempo impacto ',num2str( tempo_impact )] )