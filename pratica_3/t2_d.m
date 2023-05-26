clear all
close all

x0 = 1;v0 = 0;K = 16;m = 1;
t0=0;tf=10;

w = sqrt(K/m);

fx = @(V) V; %derivada posicao
fv = @(X) -K*X/m; %derivada tempo

h= 10.^(-1:-1:-4);
Erro = nan(1,length(h));
for i=1:length(h)
    
    t = t0:h(i):tf;
    N = length(t);
    x = nan(1,N);
    v = nan(1,N);
    x(1) = x0;v(1) = v0;

    for k=1:N-1
        r1v = fv( x(k) );
        r1x = fx( v(k) );

        r2v = fv( x(k) + r1x*h(i)/2 );
        r2x = fx( v(k) + r1v*h(i)/2 );

        r3v = fv( x(k) + r2x*h(i)/2 );
        r3x = fx( v(k) + r2v*h(i)/2 );

        r4v = fv( x(k) + r3x*h(i) );
        r4x = fx( v(k) + r3v*h(i) );

        v(k+1) = v(k) + ( r1v+2*r2v+2*r3v+r4v )*h(i)/6;
        x(k+1) = x(k) + ( r1x+2*r2x+2*r3x+r4x )*h(i)/6; 
    end
 
    x_exact = x0*cos(w.*t(N));
    Erro(i) = abs( x(N)- x_exact );
end
plot( log(h),log(Erro),'o')
lsline
%erro do metodo
aux = polyfit(log(h),log(Erro),1);
aux(1)