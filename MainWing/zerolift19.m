function [alpha0015,alpha0015m] = zerolift19(Re,codemax,U,nu)

%Re:ƒŒƒCƒmƒ‹ƒY”
%alpha00:ƒ[ƒŒ}Šp

x = [200000 250000 300000 350000 400000 450000 500000 550000 600000 650000 700000];
z = [-6.26 -7.14 -7.508 -7.842 -7.966 -8.029 -8.180 -8.242 -8.284 -8.314 -8.534];

alpha0015 = interp1(x,z,Re,'spline');

alpha0015m = interp1(x,z,U*codemax/nu,'spline');