clear all
close all

L = 1;
h = 0.025;
N_max_iter = 1E4;
tol = 1E-7;
delta_x = h;
delta_y = h;

y = -L:delta_y:L;
x = -L:delta_x:L;

%[X,Y] = meshgrid(x,y);
N = length(x);
V_old = zeros(N,N);
ind_min = L/(2*h)+1;
ind_max = 3*L/(2*h)+1;

V_old( ind_min:ind_max , ind_min:ind_max ) = 1;
figure(1)
mesh(x,y,V_old)

cond_fronteira = V_old == 0;

V_new = V_old;

for iter = 1:N_max_iter
    for i = 2:(length(x)-1)
        for j = 2:(length(y)-1)
            if true(cond_fronteira(j,i))     
                V_new(j,i) = 1/4*( V_old(j+1,i) + V_old(j-1,i) + V_old(j,i+1) + V_old(j,i-1) );
            end        
        end
    end
    num = sqrt(sum(sum((V_new-V_old).^2)));
    den = sqrt(sum(sum(V_new.^2)));
    if (num/den) < tol
        break
    end
    V_old = V_new;
end
figure(2)
mesh(x,y,V_new)

%exercicio 8.3
[Ex,Ey]=gradient(V_new,h,h);
Ex = -Ex;
Ey = -Ey;
figure(3)
quiver(x,y,Ex,Ey)

Cap_1 = -trapz(x, Ex(:,1))*4;

E2 = Ex.^2+Ey.^2;
Cap_2 = mean(mean(E2))*4;

