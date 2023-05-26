clear all
close all

h_x = 0.05;
N = 2^10;
x(:,1) = -(N-1)/2*h_x : h_x : (N-1)/2*h_x; % x centrado em 0
zi = 0;
zf = 4;
h_z = 0.02;
z(:,1) = zi : h_z : zf;
N_z = length(z);

q = zeros(N_z,N);
%a
q(1,:) = exp(-x.^2/2);
%b
q(1,:) = sech(x);

h_k = 2*pi / (N*h_x);%centrar k
k_max = (N/2-1)*h_k;
k_min = (-N/2)*h_k;
k = k_min:h_k:k_max;

q_t0 = fftshift( fft( q(1,:) ) );%fazer a transformada do q
q_t = zeros(N_z,N);

for index_z = 1:N_z
    q_t = q_t0.*exp(-1i.*k.^2.*z(index_z)/2);
    q(index_z,:) = ifft( ifftshift(q_t) );%a inversa
end

perfil = abs(q).^2;
mesh(x,z,perfil),xlabel('x'),ylabel('z'),zlabel('intensidade |q|^2')