clear all
close all
% nao faco puta ideia
L = 3; T = 5E4; w = 1E5; alfa = 5E-8;
y0 = 0;y_end = 0; B=0;
v0 = -1;
h = 0.01;
x = 0:h:L;

guess = [-0.10 -0.01];
result = nan(1,2);

reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = abstol_1;

options = odeset( 'RelTol',reltol,'AbsTol',[abstol_1 abstol_2] );
for i = 1:1:2
    for count = 1:lenght(x)
        [t_ode45,sol] = ode45( @e_2015_2016_exer_1_func,[0 L],[y0 v0],options,alfa,T,w,x(count),L);
        result(i) = sol(1,end);%ver a deflexao maxima na ponta
    end
end
m = ( result(2)-result(1) )/( guess(2) -guess(1) );
b = result(2) - m * guess(2);
guess(3) = guess(2) + ( B - result(2) )/m;

while abs(guess(2) - guess(1)) > 1E-4
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    
    [t_ode45,sol] = ode45( @e_2015_2016_exer_1_func,[0 L],[y0 v0],options,alfa,T,w,x(1),L);
    result(2) = sol(1,end);
    
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    guess(3) = guess(2) + (B-result(2))/m;
end

