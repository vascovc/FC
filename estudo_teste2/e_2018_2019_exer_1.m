clear all
close all

t0 = 0;tf = 51.1;
N = 512;
%a
delta_t = (tf-t0)/N;
f_max = 1/(2*delta_t);

%b
t = 0:delta_t:(N-1)*delta_t;
delta_w = 2*pi/(N*delta_t);

w_max = (N/2-1)*delta_w;
w_min = (-N/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;

y = 1.2*cos(5.0*t)+cos(5.1*t)+sin(5.5*t);
Z = fft(y);
Z = fftshift(Z);

dens_espt = (delta_t*abs(Z)).^2;
plot(w,dens_espt,'.')
%c
invers = ifft(ifftshift(Z));
plot(t,abs(invers-y))
%d
plot(w,dens_espt,'.')
axis([4.8 5.8 0 700])