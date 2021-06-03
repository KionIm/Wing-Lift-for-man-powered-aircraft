function [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y)

%F:�g�͙��f�́Anpartition:�������Asp�������ʒu�Aikev:�P�u���[�̈ʒu�ԍ��AWD:�㔽�p�Amb:���̋Ǐ��d�ʁAIz:�f�ʓ񎟃��[�����g�AE:�����O��
%MO:���[�����g�Atheta;����݊p�Aty:����ݗ�
%M:�Ǐ����[�����g
%Fi:���ɐ����ȉ׏d

Mij = zeros(npartition,npartition);
Fi = zeros(1,npartition);

for i = 1:npartition
    if i <= q33
        Fi(1,i) = F(1,i) - mw(1,i);
    else
        Fi(1,i) = F(1,i) - mw(1,i)*cos(thetaC);
    end
end

for i = 1:npartition
    for j = i+1:npartition
        Mij(i,j) = Fi(1,j)*(j-i)*nd;
    end
end

M = zeros(1,npartition);

for i = 1:npartition
    
    if i >= WW(1,2)+1
        M(1,i) = sum(Mij(i,:));
        
    elseif i <= WW(1,2)
        M(1,i) = sum(Mij(i,:)) - WW(1,1)*(WW(1,2)-i)*nd; 
    end
end


%����݂̎�(M/EI)�v�Z��1/EI�̍s��reEI����
reEI = zeros(1,npartition);

for i = 2:npartition
  
  reEI(1,i) = 1/(E(1,i-1)*Iz(1,i-1));
  
end
 
reEI(1,1) = 1/(E(1,1)*Iz(1,1));
 

%����݊p�v�Z
thetaij = zeros(npartition,npartition);

for i = 1:npartition 
    for j = 1:i-1
        thetaij(i,j) = M(1,j)*reEI(1,j);
    end
end

theta = zeros(1,npartition);

for i = 1:npartition    
   if i <= q33
           
        theta(1,i) = sum(thetaij(i,:)*nd);  
   else
        
        theta(1,i) = sum(thetaij(i,:)*nd) + thetaC;
   end
end



%����ݗʌv�Z
tyij = zeros(npartition,npartition);

for i = 1:npartition 
    for j = 1:i-1
        tyij(i,j) = M(1,j)*reEI(1,j)*(i-j)*nd;
    end
end    
        


ty = zeros(1,npartition);

for i = 1:npartition
    if i <= q33
           
       ty(1,i) = sum(tyij(i,:)*nd);  
    else
        
       ty(1,i) = sum(tyij(i,:)*nd) + y(1,i);
    end
end    



        
        
        