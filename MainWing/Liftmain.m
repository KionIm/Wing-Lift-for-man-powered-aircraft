function [po,effi] = Liftmain(chord,nu,alpha0015m,...
                         npartition,original15m,alpha,nd,rho,U,alpha0015,Re,codemax,alpham,...
                         bd,bd1,WW,mw,E,q33,thetaC,y0,original15mm)
%初期揚力
[F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
%抗力
[drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

%たわみ計算
[Iz] = areainertia19(bd,bd1);
[M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y0);
Fm = F.*cos(theta);

x = zeros(1,npartition);
for i = 1:npartition
     if i <= q33
          x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
     else 
          x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
     end
end
y = ty;

for i = q33+1:npartition
        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));
end

%翼中央の循環
gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
%循環分布計算
[gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
%吹き下し計算
[X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
%循環分布再計算
[gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
%局所揚力再計算
F = rho*U*gammas.*cos(theta)*nd;
%局所誘導抗力
dD = rho*downwash.*gammas;

effi = (sum(D) - sum(dD))/sum(F);

po = (sum(D) - sum(dD))*U; 

end