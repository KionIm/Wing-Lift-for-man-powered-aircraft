function [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015)
%DAE31m:�g�͌X�΁Aalpha:�݌v�}�p�Acy:�R�[�h���And:�ŏ������P�ʁArho:�C�̖��x�AU:�@���Amb:���̋Ǐ��d��
%F:�g��

F = zeros(1,npartition);

for i = 1:npartition 
     
        F(1,i) = 1/2*original15m(1,i)*chord(1,i)*(alpha(1,i)-alpha0015(1,i))*rho*U^2*nd;
 
end    
