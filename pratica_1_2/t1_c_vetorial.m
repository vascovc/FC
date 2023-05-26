clear all
close all
%em coluna e mais indicado
v0 = 50; ang = 37; alt = 0;
g = [0 0 9.8];
h = 0.01;
t = 0:h:10;
N = length(t);
pos = nan(N,3);
vel = nan(N,3);
vel(1,1) = v0*cosd(ang); %v_x
vel(1,3) = v0*sind(ang); %v_z
pos(1,1) = 0;            %x
pos(1,3) = alt;          %z

for k=1:N-1
    if pos(k,3)<0 %z<0
        break
    end
    vel(k+1,:) = vel(k,:)+(-g).*h;
    pos(k+1,:) = pos(k,:)+vel(k,:)*h;
end
t_impacto = interp1([pos(k-1,3) pos(k,3)], [t(k-1) t(k)],0)
alcance = interp1([t(k-1) t(k)], [pos(k-1,1) pos(k,1)],t_impacto)