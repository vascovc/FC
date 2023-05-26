function [ derivadas ] = f(t,sol,B)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

derivadas = zeros(3,1);
x = sol(1);
v = sol(2);
a = sol(3);

derivadas(1) = v;
derivadas(2) = a;
derivadas(3) = -0.3*a-v-B*x+sign(x);
end
