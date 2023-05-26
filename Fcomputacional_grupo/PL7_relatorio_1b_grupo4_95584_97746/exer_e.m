%PL7
%Tiago Santos 95584
%Vasco Costa  97746
clear all
close all

L = 4;
h = 0.025;
N_max_iter = 1E6;
tol = 1E-7;

y = -L/2:h:L/2;
x = -L/2:h:L/2;

N = length(x);
z_old = zeros(N,N);

z_old(1,:) = 10;
z_old(end,:) = 12;
z_old(:,1) = 11+2/L.*x;
z_old(:,end) = 11+2/L.*x;

[~,ind] = min(abs(x-0));%metodo geral, mais util para se nao fosse x=0 a querer-se

z_old(ind,ind) = 11;%a vareta
f = zeros(N,N);
for j = 1:N
    for i = 1:N 
        if (y(i)^2+ x(j)^2) < (L^2/9)
            f(j,i) = 1;
        end
    end
end

cond_front = z_old == 0; %para a vareta lá ficar colocada

z_new = z_old;
alfa = 2-2*pi/N;
for iter = 1:N_max_iter
    for i = 2:(N-1)
        for j = 2:(N-1)
            if true(cond_front(j,i)) % para nao se alterar as condicoes fronteiras
                z_new(j,i) = (1-alfa)*z_old(j,i)+alfa/4*( z_new(j+1,i) + z_new(j-1,i) + z_new(j,i+1) + z_new(j,i-1) - h^2*f(j,i) );        
            end
        end
    end
    num = sqrt(sum(sum((z_new-z_old).^2)));
    den = sqrt(sum(sum(z_new.^2)));
    if (num/den) < tol
        break
    end
    z_old = z_new;
end

[z_x,z_y]=gradient(z_new,h,h);

figure(1)
quiver(x,y,z_x,z_y),xlabel('x'),ylabel('y'),title('Vara fina, gradiente, sobre-relaxação sucessiva')
axis equal
axis([-2 2 -2 2])

mod = (z_x.^2+z_y.^2).^(1/2);
figure(2)
plot(y,z_old(:,ind)),xlabel('y'),ylabel('z'),title('z em função de y, x=0')

figure(3)
%ind dá o valor em que x=0
plot(y,mod(ind,:))% mod(:,ind) da a coluna de x ind
xlabel('y'),ylabel('modulo do gradiente'),title('modulo do gradiente de z em função de y, x=0')
