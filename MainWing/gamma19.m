function[gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U)

%gamma:èzä¬ï™ïz

gamma = zeros(1,npartition);

for i = 1:npartition 
    
        gamma(1,i) = 1/2*original15m(1,i)*chord(1,i)*(alpha(1,i)-alpha0015(1,i))*U;
    
end   

