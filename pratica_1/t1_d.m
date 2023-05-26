clear all
clc
close all

m = 0.150; g = 9.8; v0 = 0;z0 =3*2;
%velocidade
t=0:0.2:1.5;
N = length(t);v=zeros(1,N);
v(1) = v0;

for k=1:(N-1) %calcula o ponto a frente por isso so vamos ate N-1
    v(k+1) = v(k)+ (-g)*0.2;
end
%alinea e) para baixo

v_anali = velocity(g,t,v0);
plot(t,v,'x--',t,v_anali,'o-')

function v = velocity(g,t,v0)
    v = -g.*t+v0;
end