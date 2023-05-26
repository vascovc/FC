clear all
close all

%a
omega = 7.292E-5; fi = 48.865; g=9.8;L = 67;
x0 = 2;y0=0;v0=0;
h = 0.01;tf = 1.1441e+05;
t = 0:h:tf;
N = length(t);

x = nan(1,N);y=nan(1,N); v_x = nan(1,N);v_y = nan(1,N);
x(1)=x0;y(1)=y0;
v_x(1) = v0;v_y(1) = v0;
for k=1:N-1
    v_x(k+1) = v_x(k) + ( 2*omega*sind(fi)*v_y(k)-g/L*x(k) ) * h;
    v_y(k+1) = v_y(k) + (-2*omega*sind(fi)*v_x(k)-g/L*y(k) ) * h;
    
    x(k+1) = x(k) + v_x(k+1)*h;
    y(k+1) = y(k) + v_y(k+1)*h;
end
plot(x,y,'.')

%b
i=1;
for k=2:N-1
    if x(k-1)<x(k) && x(k)> x(k+1)
        teta(i) = atan( y(k)/x(k) );
        tempo(i) = t(k);
        i = i+1;
    end
end
plot(tempo,teta)

aux = polyfit(tempo,teta,1);
periodo = 2*pi/abs(aux(1))
