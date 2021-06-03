function [mb] = mbeam(m0beam,npartition,q1,q2,q3,q4,q5,q6)

%npartition:分割数、robeam:桁密度、m0beam:桁重量[1,4]
%mb:桁の局所重量
%q11,q22,q33,q44:桁分割における１、２、３、４番桁までの分割数

mmbeam = zeros(1,npartition);
mb = zeros(1,2*npartition+1);

for i = 1:q1
  
  mmbeam(1,i) = m0beam(1,1)/(2*q1);
  
end 

for i = q1+1:q2
  
  mmbeam(1,i) = m0beam(1,2)/(q2-q1-1);
  
end

for i = q2+1:q3
  
  mmbeam(1,i) = m0beam(1,3)/(q3-q2-1);
  
end

for i = q3+1:q4
    
  
  mmbeam(1,i) = m0beam(1,4)/(q4-q3-1);
  
end  

for i = q4+1:q5
    
  
  mmbeam(1,i) = m0beam(1,5)/(q5-q4-1);
  
end  

for i = q5+1:q6
    
  
  mmbeam(1,i) = m0beam(1,6)/(q6-q4-1);
  
end  
mb = mmbeam;

