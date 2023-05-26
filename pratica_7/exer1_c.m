clear all
close all

delta_t = 0.05;
N = 2^15;
t = 0:delta_t:(N-1)*delta_t;
delta_w = 2*pi/(N*delta_t); %tem que se diminuir o delta_w

w_max = (N/2-1)*delta_w;
w_min = (-N/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;

y = sin(10*t)+sin(10.05*t);
%plot(t,y),title('original')

Z = fft(y);
Z = fftshift(Z);

dens_espt = (delta_t*abs(Z)).^2;

figure(2)
plot(w,dens_espt,'.')

n_quist = pi/delta_t;

%alinea d
clear all

delta_t = 0.1;
N = 2^10;
t = 0:delta_t:(N-1)*delta_t;
delta_w = 2*pi/(N*delta_t); %tem que se diminuir o delta_w

w_max = (N/2-1)*delta_w;
w_min = (-N/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;

y = exp(-10*1i*t)+exp(20*1i*t);
%plot(t,y),title('original')

Z = fft(y);
Z = fftshift(Z);

dens_espt = (delta_t*abs(Z)).^2;

figure(2)
plot(w,dens_espt,'.')

n_quist = pi/delta_t;
