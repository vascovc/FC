clear all
close all

delta_t = 0.05;%muda se isto por causa do n_quist que dava 31, nao ia ver o 40
N = 2^12;
t = 0:delta_t:(N-1)*delta_t;
delta_w = 2*pi/(N*delta_t);

w_max = (N/2-1)*delta_w;
w_min = (-N/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;

y = sin(10*t)+sin(40*t);
%plot(t,y),title('original')
Z = fft(y);
Z = fftshift(Z);

dens_espt = (delta_t*abs(Z)).^2;

figure(2)
plot(w,dens_espt,'.')