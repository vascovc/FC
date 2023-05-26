clear all
close all

M=1.5;K=2;
t0 = 0;tf = 50;
x0 = 1.9;
v0 = 0;
A_final = -1.5;%vai ser o valor do B

alfa = -0.2;

guess = [alfa -0.15];
result = nan(1,2);

%ode45
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = abstol_1;
options = odeset( 'RelTol',reltol,'AbsTol',[abstol_1 abstol_2] );

for i=1:1:2
    [t_ode45,sol] = ode45( @f,[t0 tf],[x0 v0],options,M,K,guess(i));

    for k=1:(length(t_ode45)-1)%para encontrar o minimo
        if sol(k+1,1) > sol(k,1)
            amplitude = (sol(k,1) + sol(k+1,1) )/2;
            result(i) = amplitude;% conseguimos assim ver a estimativa para a amplitude minima
            break
        end
    end % tambem era possivel fazer uso do min neste caso
end
m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
b = result(2)-m*guess(2);
B = A_final;
guess(3) = guess(2) + (B-result(2))/m;

while abs(guess(2) - guess(1)) > 1E-3
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    %
    [t_ode45,sol] = ode45( @f,[t0 tf],[x0 v0],options,M,K,guess(2));
    
    for k=1:(length(t_ode45)-1)
        if sol(k+1,1) > sol(k,1)
            amplitude = (sol(k,1) + sol(k+1,1) )/2;
            result(2) = amplitude;
            break
        end
    end
    %
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    guess(3) = guess(2) + (B-result(2))/m;
end
disp(['alfa ',num2str(guess(3))])
plot(t_ode45,sol(:,1))