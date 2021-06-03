function [mb] = mbeam19(m0beam,npartition,q11,q22,q33,q44,q55,nd)

%npartition:�������Arobeam:�����x�Am0beam:���d��[1,4]
%mb:���̋Ǐ��d��(kg/m*m/s^2)=(N/m)
%q11,q22,q33,q44:�������ɂ�����P�A�Q�A�R�A�S�Ԍ��܂ł̕�����

mb = zeros(1,npartition);

for i = 1:q11
  
  mb(1,i) = m0beam(1,1)*nd*10^3;
end 

for i = q11+1:q22
  
  mb(1,i) = m0beam(1,2)*nd*10^3;
  
end

for i = q22+1:q33
  
  mb(1,i) = m0beam(1,3)*nd*10^3;
  
end

for i = q33+1:q44
  
  mb(1,i) = m0beam(1,4)*nd*10^3;

end

for i = q44+1:q55
  
  mb(1,i) = m0beam(1,5)*nd*10^3;
  
end 

for i = q55;1: npartition
    
   mb(1,i) = m0beam(1,3)*nd*10^3;
  
end

