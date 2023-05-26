clear all
close all

load 'data.txt'
Fs = 1000; %frequencia de amostragem
delta_T = 1/Fs; % o periodo, ou seja 1/freq
N = length(data); % o numero de elementos
t = (0:N-1)*delta_T; % vetor de tempo
%plot(t,data),xlabel('tempo(s)'),ylabel('sinal')

f = Fs*(0:(N/2))/N;%

f_sim = nan(N,1);
f_sim(1:(N/2)+1) = -f(end:-1:1);
f_sim( (N/2):N) = f(1:1:end);

Z = fft(data);
Z = fftshift(Z);

plot(f_sim,abs(Z)),xlabel('frequencia (Hz)'),ylabel('valor da transformada')
figure(2)
dens_esp = (delta_T*abs(Z)).^2;

plot(f_sim,dens_esp)% o pedido
xlabel('frequencia (Hz)')
ylabel('densidade espetral')

med = mean(dens_esp(dens_esp<0.04));%menor que 0.04 e ruido

%estudo
index = find(dens_esp<0.02);
dens_esp(index) = 0;
figure(3)
plot(f_sim,dens_esp)
figure(4)
Z(index) = 0;
plot(t,data,t,ifft(ifftshift(Z)))
