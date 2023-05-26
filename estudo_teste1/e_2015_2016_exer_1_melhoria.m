clear all
close all

w0 = 1; q = 0.5; wd = 2/3; teta0 = 0.2; v_teta0 = 0;
tf = 100;
h = 0.01;
t = 0:h:tf;
N = length(t);
F0 = [0 0.1 1.2];

for ind_F = 1: length(F0)
    teta = nan(1,N);
    v_teta = nan(1,N);
    teta(1) = teta0;
    v_teta(1) = v_teta0;
    F = F0(ind_F);
    
    f_teta = @(V) V;
    f_v_teta = @(tempo,T,V) -w0*sin(T)-q*V+F*sin(wd*tempo);

    for k=1:N-1
        r1_t = f_teta(v_teta(k) );
        r1_v = f_v_teta( t(k),teta(k),v_teta(k) );

        r2_t = f_teta(v_teta(k) + r1_v*h/2 );
        r2_v = f_v_teta( t(k)+h/2,teta(k) + r1_t*h/2,v_teta(k) + r1_v*h/2 );

        r3_t = f_teta( v_teta(k) + r2_v*h/2 );
        r3_v = f_v_teta(t(k)+h/2, teta(k) + r2_t*h/2,v_teta(k) + r2_v*h/2 );

        r4_t = f_teta( v_teta(k) + r3_v*h );
        r4_v = f_v_teta(t(k)+h, teta(k) + r3_t*h,v_teta(k) + r3_v*h );

        teta(k+1) = teta(k) + h/6*(r1_t + 2*r2_t + 2*r3_t + r4_t);
        v_teta(k+1) = v_teta(k) + h/6*(r1_v + 2*r2_v + 2*r3_v + r4_v);
    end
    figure(1)
    plot(t,teta)
    hold on
    figure(2)
    plot(teta,v_teta)
    hold on
end
