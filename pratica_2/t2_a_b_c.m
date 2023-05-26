clear all

h = 0.001;
t(:,1) = 0:h:1;
N = length(t);
x = nan(N,1);
y = nan(N,1);
ang = nan(N,1);
v_x = nan(N,1);
v_y = nan(N,1);
r = nan(N,1);

x(1) = 0.47;
y(1) = 0;
r(1) = norm([x(1) y(1)]);
ang(1) = 0;
v_x(1) = 0;
v_y(1) = 8.2;
%Euler-Cromer
for k = 1:(N-1)
    v_x(k+1) = v_x(k)-4*pi^2 * x(k) / (r(k)^3)*h;
    v_y(k+1) = v_y(k)-4*pi^2 * y(k) / (r(k)^3)*h;
    x(k+1) = x(k) + v_x(k+1)*h;
    y(k+1) = y(k) + v_y(k+1)*h;
    
    r(k+1) = norm([x(k+1),y(k+1)]);
    %norm calcula a norma, assim tem se o raio
    ang(k+1) = mod( atan2(y(k+1),x(k+1)),2*pi);
end

figure(1)
axis([-0.5 0.5 -0.5 0.5])
set(gca,'PlotBoxAspectRatio',[1 1 1])
plot(x,y,'k.-')

%%%%b%%%

for k = 1:N-1
    if ang(k+1)<ang(k)%quando salta, porque so vai de 0 a 2pi
        break
    end
end

disp(['periodo: ',num2str(t(k)),'anos'])

%%%%%c%%%%

index_t_meio_periodo = floor(k/2);
area = nan(index_t_meio_periodo-1,1);

for k = 1:(index_t_meio_periodo-1)
    area(k) = ( (r(k+1) + r(k) )/2)^2*(ang(k+1)-ang(k))/2;
end

figure(2)
plot(t(1:index_t_meio_periodo-1),area)
