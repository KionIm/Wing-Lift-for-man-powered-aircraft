function [x,cy,q1,q2,q3,q4,q11,q22,q33,q44,nd] = chord192(t,r,sp,cd0,npartition)

%t:�e�[�p�[�_�ʒu[1,4]�Ar:�e�[�p�[��[1,4]�Asp:���X�p��[1,4]
%,cd0:�ő�R�[�h���Anpartition:�������And:�����
%x:��Έʒu�Acy:�R�[�h���Aq1,2,3,4:�e�[�p�[�_�̔ԍ��Aq11,q22,q33,q44:�������ʒu�̔ԍ�

nd = sp(1,4)/npartition;

x= zeros(1,npartition);
cy = zeros(1,npartition);

%�ԍ�����U��

for i = 1:npartition
    
    x(1,i) = i*nd;

end

for i = 1:npartition
    
    if x(1,i) <= t(1,1)
        
        q1 = i;
        
    elseif x(1,i) <= t(1,2)
        
        q2 = i;
    
    elseif x(1,i) <= t(1,3)
        
        q3 = i;
        
    else 
        
        q4 = i;
        
    end    
end

for i = 1:npartition
    
    if x(1,i) <= sp(1,1)
        
        q11 = i;
        
    elseif x(1,i) <= sp(1,2)
        
        q22 = i;
    
    elseif x(1,i) <= sp(1,3)
        
        q33 = i;
        
    else 
        
        q44 = i;
        
    end    
end
    

%�R�[�h�����v�Z
for i = 1:npartition
    
    if x(1,i) <= t(1,1)  
        
        cy(1,i) = cd0-i*nd*(((1-r(1,1))*cd0)/sp(1,1));
        
    elseif x(1,i) <= t(1,2)
        
        cy(1,i) = cy(1,q1)-(i-q1)*nd*((1-r(1,2))*cy(1,q1)/(sp(1,2)-sp(1,1)));
        
    elseif x(1,i) <= t(1,3)
        
        cy(1,i) = cy(1,q2)-(i-q2)*nd*((1-r(1,3))*cy(1,q2)/(sp(1,3)-sp(1,2)));
        
    else
        
        cy(1,i) = cy(1,q3)-(i-q3)*nd*((1-r(1,4))*cy(1,q3)/(sp(1,4)-sp(1,3)));
    
    end

end

cyy = [fliplr(cy) cd0 cy];
cy = cyy;


for i = 1:2*npartition+1
    
    x(1,i) = nd*(i-1);

end