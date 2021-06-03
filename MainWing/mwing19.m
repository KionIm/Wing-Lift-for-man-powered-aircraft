function [mw] = mwing19(m0wing,npartition,q1,q2,q3,q4,q5,q6,nd)
%mw:óÉã«èäèdó (kg/m*m/s^2)=(N)

mw = zeros(1,npartition);

for i = 1:q1
  
  mw(1,i) = m0wing(1,1)*nd;
end 

for i = q1+1:q2
  
  mw(1,i) = m0wing(1,2)*nd;
  
end

for i = q2+1:q3
  
  mw(1,i) = m0wing(1,3)*nd;
  
end

for i = q3+1:q4
  
  mw(1,i) = m0wing(1,4)*nd;

end

for i = q4+1:q5
  
  mw(1,i) = m0wing(1,5)*nd;
  
end 

for i = q5+1:q6
    
  mw(1,i) = m0wing(1,6)*nd;
  
end

