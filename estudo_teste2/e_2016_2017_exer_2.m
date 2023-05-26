clear all
close all

delta_w = 0.1;

freq_nquist = nan(1,10);
for N = 5:1:15
    val = 2^N;
    freq_nquist(N-4) = val*delta_w/2;        
end
plot( [5:15] ,log(freq_nquist),'.')

[~,I] = min(abs((freq_nquist-100)));
p = I+4;

N=2^p;
delta_t = 2*pi/(N*delta_w);

t = 0:delta_t:(N-1)*delta_t;
% sem shift
w_max = (N-1)*delta_w;
w_min = 0;
w(:,1) = w_min:delta_w:w_max;

% com shift
w_max = (N/2-1)*delta_w;
w_min = (-N/2)*delta_w;
w(:,1) = w_min:delta_w:w_max;

figure(3)
subplot(1,2,1)
y = cos(50*t)-2*sin(90.1*t+pi/3);

plot(t,y)
%b
y_t = fftshift( fft(y) );
y_t_dens = (abs(y_t)*delta_t).^2;
figure(2)
plot(w,y_t_dens,'.')

ax = ifft( ifftshift(y_t_dens) );
figure(3)
subplot(1,2,2)
plot(t,y-ax,'.')

