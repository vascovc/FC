clear all
close all

h = 0.001;
r = 0:h:50;
N = length(r);

n = 3;l=1;

y = zeros(1,N);
y(end) = 0;
y(end-1) = h/1000;%valor pequeno diferente de 0

g = zeros(1,N);
guess = [-0.06 -0.07];
result = nan(1,2);

for i = 1:1:2
    g = 2*(guess(i) + 1./r - l*(l+1)./(2*r.^2));
    for k = N-1:-1:3
        y(k-1) = (1+h^2/12*g(k-1))^-1*(2*(1-5*h^2/12*g(k))*y(k)-(1+h^2/12*g(k+1))*y(k+1));
    end
    y(1) = interp1(r(2:5),y(2:5),0,'spline');
    result(i) = y(1);
end
m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
b = result(2)-m*guess(2);
B = 0;
guess(3) = guess(2) + (B-result(2))/m;

while abs( guess(2)-guess(1) ) > 1E-12
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    
    g = 2*(guess(2)+1./r-l*(l+1)./(2*r.^2));
    for k = N-1:-1:3
        y(k-1) = (1+h^2/12*g(k-1))^-1*(2*(1-5*h^2/12*g(k))*y(k)-(1+h^2/12*g(k+1))*y(k+1));
    end
    y(1) = interp1(r(2:5),y(2:5),0,'spline');
    result(2) = y(1);
    
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    %b = result(2)-m*guess(2); % nao e necessario
    guess(3) = guess(2) + (B-result(2))/m;
    %seria de recomendar usar um valor limite para ver a covergencia
end

disp(['valor pratico obtido ',num2str(guess(3))])
disp(['valor teorico ',num2str( -1/2*n^-2 )])

R = y(2:end)./r(2:end);
C = trapz( r(2:end),R.^2 );
R_norm = R./C^(1/2);
s = trapz(r(2:end),R_norm.^2); % tem que dar 1
plot(r(2:end),R_norm)
