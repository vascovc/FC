clear all
close all

%% a

a = 2; b =0.74; c = 0.5; x0 = 0.8; y0 = 0.3;
h = 0.01;
tf = 100;
t = 0:h:tf;
N = length(t);
x = nan(1,N); y = nan(1,N);
x(1) = x0; y(1) = y0;

fx = @(X,Y) X*(1-X)-a*X*Y/(X+Y);
fy = @(X,Y) b*X*Y/(X+Y)-c*Y;

for k=1:N-1
    r1x = fx( x(k),y(k) );
    r1y = fy( x(k),y(k) );
    
    r2x = fx( x(k) + r1x*h/2,y(k) + r1y*h/2 );
    r2y = fy( x(k) + r1x*h/2,y(k) + r1y*h/2 );
    
    r3x = fx( x(k) + r2x*h/2,y(k) + r2y*h/2 );
    r3y = fy( x(k) + r2x*h/2,y(k) + r2y*h/2 );
    
    r4x = fx( x(k) + r3x*h,y(k) + r3y*h );
    r4y = fy( x(k) + r3x*h,y(k) + r3y*h );
    
    x(k+1) = x(k) + h/6*(r1x + 2*r2x + 2*r3x + r4x);
    y(k+1) = y(k) + h/6*(r1y + 2*r2y + 2*r3y + r4y);
end

plot(t,x)
hold on
%% b
a = 2; b =0.74; c = 0.5; x0 = 0.8; y0 = 0.3;
h = 0.01;
t0 = 0;
tf = 100;
x(1) = x0; y(1) = y0;

reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = abstol_1;

options = odeset( 'RelTol',reltol,'AbsTol',[abstol_1 abstol_2] );
[t_ode45,sol] = ode45( @e_2017_2018_exer_2_recurso_function,[t0 tf],[x0 y0],options,a,b,c);

plot(t_ode45,sol(:,1))

