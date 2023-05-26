clear all
close all

alfa = 40.6;
alt  = 200;
W    = 7.292*10^-5;
pos0 = [0 0 alt];
v0   = [0 0 0];
g0   = [0 0 -9.8];
r    = [0 0 6.37*10^6]; %alinhado com o z
w    = [-W*cosd(alfa) 0 W*sind(alfa)];
g    = g0 - cross(w, cross(w,r) );
%cross faz o produto externo

t0 = 0;tf = 10;
h = 0.01;
t = t0:h:tf;
N = length(t);
v = nan(N,3);
pos = nan(N,3);
v(1,:) = v0;
pos(1,:) = pos0;

for k=1:N-1
    if pos(k,3)<0
        break
    end
    v(k+1,:) = v(k,:)+g.*h;
    pos(k+1,:) = pos(k,:)+v(k,:).*h;
end

t_impacto = interp1( [pos(k-1,3) pos(k,3)], [t(k-1) t(k)],0)
x_imp     = interp1( [t(k-1) t(k)]        ,[pos(k-1,1) pos(k,1)],t_impacto)
y_imp     = interp1( [t(k-1) t(k)]        ,[pos(k-1,2) pos(k,2)],t_impacto)