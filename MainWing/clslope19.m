function [original15m,original15mm] = clslope19(Re,codemax,U,nu)

%original15—g—ÍŒXŽÎ
x = [200000,250000,300000,350000,400000,450000,500000,550000,600000,650000,700000];
z = [0.1112 0.1049 0.1025 0.10 0.0992 0.0989 0.0977 0.0973 0.097 0.0968 0.0951];
original15m = interp1(x,z,Re,'spline');

original15mm = interp1(x,z,U*codemax/nu,'spline');
