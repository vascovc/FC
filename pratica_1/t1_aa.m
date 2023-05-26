clear all
close all
clc

m = 0.150; g = 9.8; v0 = 0;z0 =3*2;

t=0:0.2:1.5;
v = velocity(g,t,v0);
z = position(g,t,v0,z0);
plot(t,v)
figure(2)
plot(t,z)

function v = velocity(g,t,v0)
    v = -g.*t+v0;
end
function z = position(g,t,v0,z0)
    z = -0.5*g.*t.^2+v0.*t+z0;
end