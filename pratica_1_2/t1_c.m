clear all
close all

v0 = 50; ang = 37; alt = 55;
g = 9.8;
h = 0.01;
t=0:h:10;

N = length(t);
vx = nan(1,N);vz = nan(1,N);
x = nan(1,N);z = nan(1,N);

z(1) = alt;x(1)=0;
vx(1) = v0*cosd(ang);
vz(1) = v0*sind(ang);

for k=1:N-1
    if z(k)<0
        break
    end
    vx(k+1) = vx(k); %segundo a horizontal nao ha forcas
    vz(k+1) = vz(k) + (-g)*h;
    x(k+1) = x(k) + vx(k)*h;
    z(k+1) = z(k) + vz(k)*h;
end

%interp mete-se o x, y, o valor que se poe no x para ter y
t_impacto = interp1([z(k-1) z(k)], [t(k-1) t(k)],0)
%o k atual e o primeiro negativo k-1 e o ultimo positivo
alcance = interp1([t(k-1) t(k)], [x(k-1) x(k)],t_impacto)