clear all
clc
close all

m = 0.150; g = 9.8; v0 = 0;z0 =3*2;h=0.2;

t=0:h:1.5;
N = length(t);v=zeros(1,N);z=zeros(1,N);
v(1) = v0;z(1)=z0;

for k=1:(N-1)
    v(k+1) = v(k)+ (-g)*h;
    z(k+1) = z(k)+ v(k)*h;
end

z_analit = position(g,t,v0,z0);
plot(t,z_analit,t,z);

function z = position(g,t,v0,z0)
    z = -0.5*g.*t.^2+v0.*t+z0;
end