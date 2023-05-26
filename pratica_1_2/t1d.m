clear all
close all
% ESTA COM VETORES LINHA E E MELHOR EM COLUNA
v0 = 50; ang = 37; alt = 0;
g = [0 ;0 ;9.8];
h = 0.01;
t = 0:h:10;
N = length(t);
pos = nan(3,N);
vel = nan(3,N);
vel(1,1) = v0*cosd(ang);
vel(3,1) = v0*sind(ang);
pos(1,1) = 0;
pos(3,1) = alt;
for k=1:N-1
    if pos(3,k)<0
        break
    end
    vel(:,k+1) = vel(:,k)+(-g).*h;
    pos(:,k+1) = pos(:,k)+vel(:,k)*h;
end
t_impacto = interp1([pos(3,k-1) pos(3,k)], [t(k-1) t(k)],0)
%no x fica o tempo e no y a posição, sabemos o tempo queremos a posicao
alcance = interp1([t(k-1) t(k)], [pos(1,k-1) pos(1,k)],t_impacto)

alcance_exato = v0^2/norm(g)*sind(2*ang)
erro = abs( (alcance-alcance_exato) /alcance_exato)