clear all
close all

y = csvread('som2.csv');
Fs = 22.05E3; %frequencia de amostragem
delta_t = 1/Fs;
N = length(y); % o numero de elementos
t = (0:N-1)*delta_t; % vetor de tempo

f = Fs*(0:(N/2))/N;%

f_sin = nan(N,1);
f_sin(1:(N/2)+1) = -f(end:-1:1);
f_sin( (N/2):N) = f(1:1:end);

[~,ind,val] = find(f_sin <= 1500);
Z = fft(y);
Z = fftshift(Z);