function [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015)
%DAE31m:揚力傾斜、alpha:設計迎角、cy:コード長、nd:最小分割単位、rho:気体密度、U:機速、mb:桁の局所重量
%F:揚力

F = zeros(1,npartition);

for i = 1:npartition 
     
        F(1,i) = 1/2*original15m(1,i)*chord(1,i)*(alpha(1,i)-alpha0015(1,i))*rho*U^2*nd;
 
end    
