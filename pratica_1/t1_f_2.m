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
    
    if z(k) < 0
        break
    end
end
plot(t,z);

% z(end-1:end)
% t(end-1:end)

tempo = interp1( z(end-1:end) ,t(end-1:end),0)
z = z(1:k-1);
v = v(1:k-1);
t = t(1:k-1);
