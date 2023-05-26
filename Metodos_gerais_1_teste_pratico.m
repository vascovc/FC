%% Instantes de interpolação
tempo = interp1( z(end-1:end) ,t(end-1:end),0);
% primeiro : dois pontos y (z)
% segundo  : dois pontos x (t)
% terceiro : valor do y que se pretende saber o x (y=0)
% devolvendo o valor de x   (o tempo em que y=0)

%% Forma de arranjar a reta
lsline
%desenha no grafico a reta de melhor ajuste
aux = polyfit(log(h),log(Erro),1);
% valores de x, valores de y, ordem do polinomio
% aux(1) declive
% aux(2) ordenada origem b

%% ordem de erro do metodo
aux = polyfit(log(h),log(Erro),1)
%em aux(1) fica o declive que indica a ordem

%% Metodo de Euler
t = t0:h:tf;
N = length(t);
for i=1:N-1
    v(k+1) = v(k) + (-g)*0.2;
    %y(k+1) = y(k) + f( t(k),y(k) ) * h
    %f representa a derivada
end
%erro de ordem h^1

%% Metodo de Euler-Cromer
v(k+1) = v(k) + f( t(k),y(k),v(k) ) * h;
y(k+1) = y(k) + v(k+1) * h;
%muda no y que fica com o v(k+1)
%erro de ordem h^???

%% Metodo de Euler - Implicito
y(k+1) = y(k) + f( t(k+1),y(k+1) ) * h;
%NECESSARIO USAR O LINSOLVE
%erro de ordem h

%% Metodo de Crank-Nicolson
y(k+1) = y(k) + ( f( t(k),y(k)) + f( t(k+1),y(k+1) )) *h/2;
%NECESSARIO USAR O LINSOLVE
%erro de ordem h^2
%% LINSOLVE
AZ=b;
A = [coeficiente_Y(k+1) coeficiente_V(k+1);
     coeficiente_Y(k+1) coeficiente_V(k+1)];
b = [tudo_linha(1);
     tudo_linha(2)];
Z = linsolve(A,b);
Z(1) = Y(k+1);
Z(2) = V(k+1);
%atualizar o b com os novos valores

%% Metodo de Runge-Kutta
fx = @(V) V; %derivada posicao
fv = @(X) -K*X/m; %derivada tempo
for k=1:N-1
    r1v = fv( x(k) );
    r1x = fx( v(k) );
    
    r2v = fv( x(k) + r1x*h/2);
    r2x = fx( v(k) + r1v*h/2);
    
    r3v = fv( x(k) + r2x*h/2);
    r3x = fx( v(k) + r2v*h/2);
    
    r4v = fv( x(k) + r3x*h);
    r4x = fx( v(k) + r3v*h);
    
    v(k+1) = v(k) + ( r1v + 2*r2v + 2*r3v + r4v )*h/6;
    x(k+1) = x(k) + ( r1x + 2*r2x + 2*r3x + r4x )*h/6; 
end
%erro de ordem do numero do metodo

%% ODE45
    % *funcao em ficheiro separado*
        function derivadas = f(t,solucao)
            derivadas = zeros(2,1);
            % Esta linha é necessária e faz do vetor de saída um vetor coluna.
            K=16;
            m=1;
            % O vetor solucao tem os valores de x e v para o tempo
            derivadas(1) = solucao(2);
            %derivada de x e a velocidade
            
            derivadas(2) = -K/m * solucao(1);
            %a derivada da velocidade e -k sobre m vezes x
        end
    % *codigo principal*
reltol = 3E-14;
%erro relativo

abstol_1 = 1E-13;
abstol_2 = abstol_1;
%intervalo do erro absoluto

options = odeset( 'RelTol',reltol,'AbsTol',[abstol_1 abstol_2] );
[t_ode45,sol] = ode45( @f,[t0 tf],[x0 v0],options);
% ( @func, intervalo do tempo, [x v], options, parametros extra para a func)
%devolve: [ tempo , sol(1) = x    sol(2) = v ]

