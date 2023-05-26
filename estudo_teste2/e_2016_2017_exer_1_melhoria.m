clear all
close all

L = 2;
C = 0;
val_analit = pi^2/(2*L)^2;
guess = [val_analit*0.9 val_analit*1.1];

h = 0.05;
x = 0:h:L;
N = length(x);
fi = nan(1,N);
deriv = nan(1,N);

fi(1) = 0;
deriv(1) = 4;

result = nan(1,2);

for i=1:1:2
    for k = 1:N-1
        deriv(k+1) = (guess(i)*fi(k)-C*x(k)*fi(k))/(-2);
        fi(k+1) = fi(k) + h*deriv(k+1);
    end
    result(i) = fi(end);
end
m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
b = result(2)-m*guess(2);
B = 0;
guess(3) = guess(2) + (B-result(2))/m;

while abs(guess(2) - guess(1)) > 1E-3
    guess(1) = guess(2);
    result(1) = result(2);
    guess(2) = guess(3);
    result(2) = result(1);
    for k = 1:N-1
        deriv(k+1) = (guess(i)*fi(k)-C*x(k)*fi(k))/(-2);
        fi(k+1) = fi(k) + h*deriv(k+1);
    end
    result(2) = fi(end);
    m = ( result(2)-result(1) )/( guess(2)-guess(1) ); 
    guess(3) = guess(2) + (B-result(2))/m;
end