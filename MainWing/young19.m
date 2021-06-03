function [E] = young19(npartition,E0,q11,q22,q33,q44,q55,q66,q77)

%npartition:�������AE0:���̃����O��[1,4]�Aq11,22,33:�������ʒu�ԍ�
%E:�����O��

E = zeros(1,npartition);

for i = 1:npartition
  
  if i <= q11
    
    E(1,i) = E0(1,1);
  
  elseif i <= q22
    
    E(1,i) = E0(1,2);
  
  elseif i <= q33
    
    E(1,i) = E0(1,3);
    
  elseif i <= q44
      
    E(1,i) = E0(1,4);
    
  elseif i <= q55
      
    E(1,i) = E0(1,5);
    
  elseif i <= q66
    
    E(1,i) = E0(1,6);
    
  else i <= q77
      
    E(1,i) = E0(1,7);
    
  end
  
end



