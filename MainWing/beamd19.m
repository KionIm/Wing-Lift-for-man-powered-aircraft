function [bd,bd1,br,br1] = beamd19(beamin,beamsp,npartition,q11,q22,q33,q44,q55,q66,q77,nd)

%beam:ͺ[a[1,8]Abeamsp:Xp[1,6]Anpartition:ͺ
%bd:OaAbd1:ΰaAbr:OΌaAbr1:ΰΌaAbx:ξͺΚu
%bt:e[p[

bd = zeros(1,npartition);
bd1 = zeros(1,npartition);
bx = zeros(1,npartition);
bt = zeros(1,7);

%e[p[
for i = 1:7
  
  if i == 1
   
    bt(1,i) = (beamin(1,2*i)-beamin(1,2*i-1))/beamsp(1,i);
    
  else
    
    bt(1,i) = (beamin(1,2*i)-beamin(1,2*i-1))/(beamsp(1,i)-beamsp(1,i-1));
    
  end
  
end

%ͺΚu
for i = 1:npartition
  
  bx(1,i) = i*nd;

end

%ΰa
for i = 1:q11
  
  bd1(1,i) = beamin(1,1)-bt(1,1)*bx(1,i);
  
end

for i = q11+1:q22
  
  bd1(1,i) = beamin(1,3)-bt(1,2)*(bx(1,i)-beamsp(1,1));
  
end

for i = q22+1:q33
  
  bd1(1,i) = beamin(1,5)-bt(1,3)*(bx(1,i)-beamsp(1,2));
  
end

for i = q33+1:q44
  
  bd1(1,i) = beamin(1,7)-bt(1,4)*(bx(1,i)-beamsp(1,3));
  
end

for i = q44+1:q55
  
  bd1(1,i) = beamin(1,9)-bt(1,5)*(bx(1,i)-beamsp(1,4));
  
end

for i = q55+1:q66
    
  bd1(1,i) = beamin(1,11)-bt(1,6)*(bx(1,i)-beamsp(1,5));
  
end

for i = q66+1:q77
    
  bd1(1,i) = beamin(1,13)-bt(1,7)*(bx(1,i)-beamsp(1,6));
  
end

%Oa
bd = bd1 + 0.002;


%OΌa
br = bd/2;
br1 = bd1/2;

end


    


