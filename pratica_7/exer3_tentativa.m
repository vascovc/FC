clear all
close all

h_x = 0.05;
N = 2^10;
x(:,1) = -(N-1)/2*h_x:h_x:(N-1)/2*h_x;% x centrado em 0
zi = 0;
zf = 4;
h_z = 0.02;
z(:,1) = zi:h_z:zf;
N_z = length(z);

q = zeros(N_z,N);
q(1,:) = sech(x);

h_k = 2*pi/(N*h_x);%centrar k como sempre
k_max = (N/2-1)*h_k;
k_min = (-N/2)*h_k;
k = k_min:h_k:k_max;
 
abstol=ones(1,N);
abstol=1E-9.*abstol;
options= odeset('RelTol',1e-9,'AbsTol',abstol);

q_t0 = fftshift( fft( q(1,:),N ));
q_t = q_t0.*exp(1i.*k.^2.*z./2);
[z,qtexp] = ode45(@nonlinear,[zi zf],[q_t0],options,N,k);
for j=1:numel(z)
    qz = qtexp(j,:).*exp(-1i.*k.^2.*z(j)./2);
    q(j,:) = ifft(ifftshift(qz));
end
mesh(x,z,abs(q).^2)
figure(2)
contourf(x,z,abs(q).^2)