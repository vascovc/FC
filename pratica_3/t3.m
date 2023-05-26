clear all
close all
 
h = 0.01;
x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;
t = t0:h:tf;
w = sqrt(K/m);
x(1) = x0;v(1) = v0;

reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = abstol_1;

options = odeset( 'RelTol',reltol,'AbsTol',[abstol_1 abstol_2] );
[t_ode45,sol] = ode45( @f,[t0 tf],[x0 v0],options);

subplot(1,2,1)
plot(t_ode45,sol(:,1),t,x0*cos(w.*t),'.')
subplot(1,2,2)
plot(t_ode45,sol(:,2),t,-w*x0*sin(w.*t),'.')