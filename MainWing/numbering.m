function [q1,q2,q3,q4,q5,q6,q11,q22,q33,q44,q55,q66,q77] = numbering(tc1,tc2,tc3,tc4,tc5,tc6,beamsp)

%t:�e�[�p�[�_�ʒu[1,4]�Ar:�e�[�p�[��[1,4]�Asp:���X�p��[1,4]
%npartition:�������And:�����
%q1,2,3,4:�e�[�p�[�_�̔ԍ��Aq11,q22,q33,q44:�������ʒu�̔ԍ�

q1 = round(tc1/0.025);

q2 = round(tc2/0.025);

q3 = round(tc3/0.025);

q4 = round(tc4/0.025);

q5 = round(tc5/0.025);

q6 = round(tc6/0.025);

q11 = round(beamsp(1,1)/0.025);

q22 = round(beamsp(1,2)/0.025);

q33 = round(beamsp(1,3)/0.025);

q44 = round(beamsp(1,4)/0.025);

q55 = round(beamsp(1,5)/0.025);

q66 = round(beamsp(1,6)/0.025);

q77 = round(beamsp(1,7)/0.025);


