clear all
close all

L = 1;
D = 0.02; m=5;
tf = 3;
delta_x = 0.01;
delta_t = 0.001;
x = 0:delta_x:L;N_x = length(x);
t = 0:delta_t:tf;N_t = length(t);

u = zeros(N_x,N_t);

u(:,1) = 0.5;
u((x<0.6),1) = 0;
u((x>0.8),1) = 0;

for ind_t=1:N_t-1
    for ind_x = 2:N_x-1
       u(ind_x,ind_t+1) = u(ind_x,ind_t)+delta_t*(D*( u(ind_x+1,ind_t)-2*u(ind_x,ind_t)+u(ind_x-1,ind_t))/delta_x^2+m*u(ind_x,ind_t)*(1-u(ind_x,ind_t)));
    end
end

figure(1)
mesh(t,x,u)
xlabel('t')
ylabel('x')
zlabel('u')
figure(2)
contourf(x,t,u')
colorbar