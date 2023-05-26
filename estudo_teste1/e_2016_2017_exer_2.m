clear all
close all

L = 0.1; ro = 2.70E3; stefan = 5.6703E-8;C = 0.91E3;
T0 = 2000; t_f = 1000;
h=0.01;
t=0:h:t_f;
N = length(t);

A = 6 * L^2;
m = ro * L^3;

c(1) = -stefan*A*h/(2*m*C);
c(2:3) = 0;
c(4)=1;c(5)=0;
sol = roots(c);
[valor, indice] = min(abs(sol-T(k)));
T(k+1)=sol(indice);
for k=1:N-1
    
end

%b
h = .1:E-1:.1E-9;

tf = 1000;
len_h = length(h);
erro(index_h) = nan(1,len_h);
for index_h=1:1:len_h
    T_exact = ( T0^-3 + 3*stefan*A/(m*C)*tf)^-1/3;
    %usar metodo anterior ate 1000
    ...
    erro(index_h) = T_exact - T(k);
end
plot(log(h),erro)

%c
L = 0.1; ro = 2.70E3; stefan = 5.6703E-8;C = 0.91E3;
t_f = 60*2*60;T0 = 310;
h=0.01;
t=0:h:t_f;
N = length(t);

T = nan(1,N);
T(1) = T0;
A = 6 * L^2;
m = ro * L^3;

%Tc = 283+1.0E-3*t;
f_Tc = @(tempo) 283+tempo*1.0E-3;
f_T = @(tempo,Temp) -stefan*A/(m*C)*(Temp.^4- (f_Tc(tempo)).^4) ;
for k=1:N-1
    r1 = f_T( t(k),T(k) ); 
    
    r2 = f_T( t(k)+h/2, T(k) + r1*h/2 );
    
    r3 = f_T( t(k)+h/2, T(k) + r2*h/2 );
    
    r4 = f_T( t(k)+h, T(k) + r3*h  );
    
    T(k+1) = T(k) + h/6*(r1 + 2*r2 + 2*r3 + r4);
end
plot(t,T)

%b
L = 0.1; ro = 2.70E3; stefan = 5.6703E-8;C = 0.91E3;
T0 = 2000; t_f = 1000;

A = 6 * L^2;
m = ro * L^3;
h = 10.^(-1:-1:-4);

len_h = length(h);
erro = nan(1,len_h);
f_T = @(tempo,Temp) -stefan*A/(m*C)*(Temp.^4) ;
for index_h=1:1:len_h
    t=0:h(index_h):t_f;
    N = length(t);
    T = nan(1,N);
    T(1) = T0;
    for k=1:N-1
        r1 = f_T( t(k),T(k) ); 

        r2 = f_T( t(k)+h(index_h)/2, T(k) + r1*h(index_h)/2 );

        r3 = f_T( t(k)+h(index_h)/2, T(k) + r2*h(index_h)/2 );

        r4 = f_T( t(k)+h(index_h), T(k) + r3*h(index_h)  );

        T(k+1) = T(k) + h(index_h)/6*(r1 + 2*r2 + 2*r3 + r4);
    end
    T_exact = ( T0^-3 + 3*stefan*A/(m*C)*t_f)^-1/3;
    erro(index_h) = T_exact - T(k);
end
plot(log(h),erro)