clear
clc
close

zi = 0;
zf = 4;
z = zi:0.02:zf;
dx = 0.1;
x = -51.15:dx:51.15;
N = numel(x);
dk = 2*pi/(N*dx);
k = -N/2*dk:dk:(N/2-1)*dk;

phi(1,:) = sech(x);
qt0 = abs(fftshift(fft(phi(1,:),N)));

abstol=ones(1,N);
abstol=1e-9.*abstol;
options= odeset('RelTol',1e-9,'AbsTol',abstol);

[z,qt]=ode45(@nonlinear,[zi zf],[qt0],options,N,k);

for j = 1:numel(z)
    f2 = qt(j,:).*exp(-1i/2*z(j)*k.^2);
    phi(j,:) = abs(ifftshift(ifft(f2)));
end

mesh(x,z,phi)