clear all
close all

y0_vetor = [2 -5 0.2];
v0_vetor = [7 -2 0.7];
epsilon = 1;
t0=0;tf=100;

reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = abstol_1;
for index = 1:3
    options = odeset( 'RelTol',reltol,'AbsTol',[abstol_1 abstol_2] );
    [t_ode45,sol] = ode45( @f,[t0 tf],[y0_vetor(index) v0_vetor(index)],options,epsilon);
    %para passar o epsilon tem que se meter no fim de tudo
    subplot(1,3,1)
    plot(t_ode45,sol(:,1))
    hold all
    subplot(1,3,2)
    plot(t_ode45,sol(:,2))
    hold all
    subplot(1,3,3)
    plot(sol(:,1),sol(:,2));
    hold all
end