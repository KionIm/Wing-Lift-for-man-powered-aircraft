function [mb] = mbeam19(m0beam,npartition,q11,q22,q33,q44,q55,nd)

%npartition:分割数、robeam:桁密度、m0beam:桁重量[1,4]
%mb:桁の局所重量(kg/m*m/s^2)=(N/m)
%q11,q22,q33,q44:桁分割における１、２、３、４番桁までの分割数

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

