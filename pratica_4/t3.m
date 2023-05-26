clear all
close all

y0 = 0.2;
v0 = 0.7;
epsilon = 0.1;
t0=0;tf=100;
y(1) = y0;v(1) = v0;

reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = abstol_1;

options = odeset( 'RelTol',reltol,'AbsTol',[abstol_1 abstol_2] );
[t_ode45,sol] = ode45( @f,[t0 tf],[y0 v0],options,epsilon);
%para passar o epsilon tem que se meter no fim de tudo
subplot(1,3,1)
plot(t_ode45,sol(:,1))
subplot(1,3,2)
plot(t_ode45,sol(:,2))
subplot(1,3,3)
plot(sol(:,1),sol(:,2));