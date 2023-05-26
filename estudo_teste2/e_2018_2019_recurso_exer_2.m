clear all
close all

load 'sinal.txt'
delta_t = 0.01;
N = length(data); % o numero de elementos
t = (0:N-1)*delta_T;

f = 1/delta_t*(0:(N/2))/N;

w_sin = nan(N,1);
w_sin(1:(N/2)+1) = -f(end:-1:1);
w_sin( (N/2):N) = f(1:1:end);

Z = fft(data);
Z = fftshift(Z);
plot(w_sin,abs(Z)),xlabel('frequencia (Hz)'),ylabel('valor da transformada')
figure(2)
dens_esp = (delta_T*abs(Z)).^2;

plot(f_sim,dens_esp)% o pedido
xlabel('frequencia (Hz)')
ylabel('densidade espetral')

%b
