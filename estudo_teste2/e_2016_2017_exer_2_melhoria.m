clear all
close all

delta_w = 0.1;
N = 1024;
delta_t = 2*pi/(N*delta_w);

t = 0:delta_t:(N-1)*delta_t;
y = cos(2*t)+0.5*cos(4*t)+1/6*cos(12*t);

w_max = (N/2-1)*delta_w;
w_min = (-N/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;

Z = fftshift(fft(y));
dens_espt = (delta_t*abs(Z)).^2;
plot(w,dens_espt)

%b
clear all
close all

load 'temperatura.txt'

delta_t = 1/12;
Fs = 12;

N = length(temperatura);
t = (0:N-1)*delta_t;
plot(t,temperatura)
f = Fs*(0:(N/2))/N;%

f_sim = nan(N,1);
f_sim( 1:((N/2)+1) ) = -f(end:-1:1);
f_sim( (N/2):N) = f(1:1:end);

Z = fft(temperatura);
Z = fftshift(Z);
dens_esp = (delta_t*abs(Z)).^2;
plot(f_sim,dens_esp )

index = find(dens_esp < 1);
dens_esp(index)=0;
plot(f_sim,dens_esp )

inver = ifft(ifftshift(dens_esp));
plot(t,inver,t,temperatura)