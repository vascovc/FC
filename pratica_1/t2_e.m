clear all
close all
clc

K = 2.5; m = 0.5; x0 = 0.10 ;v0 = 0;

h = [ 1E-1,5E-2,1E-2,5E-3,1E-3,5E-4,1E-4 ];
%h = 10.^(-1:-1:-4);
erro = nan(1,length(h));

for index_h = 1:length(h)
    t = 0:h(index_h):10 ;
    N = length(t);
    v = zeros(1,N); x = zeros(1,N);
    v(1) = v0;x(1) = x0;
    for k = 1:N-1 %Euler
        v(k+1) = v(k) + (-K/m)*x(k)*h(index_h);
        x(k+1) = x(k) + v(k)*h(index_h);
    end
    erro(index_h) = abs(x0-max(x));
end

plot(log10(h),log10(erro),'o')
lsline %da a reta de melhor ajuste
aux = polyfit(log10(h),log10(erro),1); %obter a reta (1)declive (2)ordenada_origem
declive = aux(1)