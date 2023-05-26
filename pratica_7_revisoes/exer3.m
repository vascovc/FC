clear all
close all

alfa = 0.2;

x_min = -20; x_max = 20;
delta_x = 0.5;
x(:,1) = x_min:delta_x:x_max;
N_x = length(x);

z_min = 0; z_max = 16;
delta_z = 0.1;
z(:,1) = z_min:delta_z:z_max;
N_z = length(z);

FI = nan(N_x,N_z);
FI(:,1) = exp(-1/2.*x.^2);
FI(1,:) = 0;
FI(N_x,:) = 0;

eta = 1i*delta_z/(4*delta_x^2);
psi = 2*alfa*delta_x^2;

A1 = diag( 1/eta + 2 + psi*( x(2:N_x-1) ).^2 );
A2 = diag( -ones(N_x-3,1),1);
A3 = diag( repmat( -1, 1, N_x-3 ),-1);
A = A1+A2+A3;

B = nan(N_x-2,1);
for n = 1:N_z-1
    for i = 2:N_x-1
        B(i-1) = FI(i-1,n) + ( 1/eta-2-psi*x(i).^2)*FI(i,n)+FI(i+1,n);
    end
    Z = linsolve(A,B);
    FI( 2:(N_x-1) ,n+1) = Z;
end

mesh(x,z,abs(FI')),xlabel('x'),ylabel('z'),colorbar
figure(2)
contourf(x,z,abs(FI')),xlabel('x'),ylabel('z'),colorbar