clear all
close all

delta_t = 0.1;
N = 2^10;
t = 0:delta_t:(N-1)*delta_t;
delta_w = 2*pi/(N*delta_t);

%% sem shift
w_max = (N-1)*delta_w;
w_min = 0;
w(:,1) = w_min:delta_w:w_max;

%% com shift
w_max = (N/2-1)*delta_w;
w_min = (-N/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;

%% resto igual
y = sin(t)+sin(10*t);
plot(t,y),title('original')
Z = fft(y);
Z = fftshift(Z);

dens_espt = (delta_t*abs(Z)).^2;

figure(2)
plot(w,dens_espt,'.')

n_quist = pi/delta_t;
%estudo
figure(3)
plot(t,y,t, ifft(ifftshift(Z)))