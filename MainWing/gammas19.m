function[gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash)

%gamma:èzä¬ï™ïz

gammas = zeros(1,npartition);

for i = 1:npartition 
    
        gammas(1,i) = 1/2*original15m(1,i)*chord(1,i)*(alpha(1,i)-alpha0015(1,i)+180/pi*downwash(1,i)/U)*U;
    
   
end   

