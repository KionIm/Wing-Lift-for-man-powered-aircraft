%%�嗃��͐݌v�v���O���� 2019�N���s���͍���H��
%{
2015,2016�N���s���͐݌v���_�c���i���񂪐��삵��MainEruba.m��2017�N���s���͐݌v�Җ؉��Y���
�����삵��DaidalousMain.m���Q�l�ɍ���Ă��܂��B
%}
%�e�[�p�[�ύX�_�ƌ������_��0.025m�̔{���ɂȂ�悤�ɂ��Ă��������B

% tc?=�e�[�p�[�ύX�_(m)
tc1=1.75;
tc2=4.75;
tc3=7.75; 
tc4=10.75;
tc5=13.75;
tc6=13.775;
%��SPAN(m)([1�ԗ�,2�ԗ�,3�ԗ�,4�ԗ�,5�ԗ�,6�ԗ�])
sp = [1.75 4.75 7.75 10.75 13.75 13.775];
%���X�p��(m)
beamsp = [1.75 4.75 7.75 10.75 13.0 13.75 13.775];
% �k=����(m)
L = sp(1,6);
WD0 = 8;%���[����݊p(�㔽�p)(deg)
thetaC = 0/180*pi;%�r���㔽�p
WD = 8/180*pi;%(rad)
W0 = 51;%�p�C���b�g�d��(kg�d)
W1 = 42;%�@�̏d��(kg�d)
W = W0+W1;%���d��(kg�d)
%�d�͉����x(m/s^2)
g = 9.797;



%----------------------------------%
% number=�����ԍ�
number = 1 : L/0.025;
%��������
npartition = round(L/0.025);
%���������
nd = sp(1,6)/npartition;
%�e�[�p�[�ύX�_�A�������_�̕����_�ԍ�
[q1,q2,q3,q4,q5,q6,q11,q22,q33,q44,q55,q66,q77] = numbering(tc1,tc2,tc3,tc4,tc5,tc6,beamsp);
%�����_�̍��W(�������̗��}�E���g��x=0,y=0)
x = zeros(1,npartition);
for i = 1:npartition
    if i <= q33 
        x(1,i) = nd*i;
    else
        x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
    end
end

y = zeros(1,npartition);
for i = 1:npartition
    if i > sp(1,3)/0.025
        y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
    end
end



%�v�Z�ϐ�----------------------------------------
% codemax=�ő嗃����(m)
codemax=1.045;
%������(�e�[�p�[�_�j
chord0 = codemax;
chord1 = 1.04;
chord2 = 0.968;
chord3 = 0.892;
chord4 = 0.786;
chord5 = 0.631;
chord6 = 0.625;

%V=�@�̑��x(m/s)
U = 7.2;

% �������͒l[���������׏d(N) ���t���ʒu(number)]
WW = [200 336];

%alpha=�ݒ�}�p(���j
alpham = 4.5;%�����}�p
alpha1 = 4.5;
alpha2 = 4.5;
alpha3 = 4.4;
alpha4 = 4.1;
alpha5 = 3.6;
alpha6 = 3.55;

%�v�Z��
N = 100;
% cr=1���s���Ƃ�chord�ω���(m)
cr = 0.001;
%---------------------------------------------------

%nu=���S���W��(m^2/s)
nu = 1.55*10^(-5);
%rho=��C���x(kg/m^3)
rho = 1.185;

% count=���s��
count = 1;


%�����a(m)([1�Ԍ���,�[,2�Ԍ���,�[,3�Ԍ���,�[,4�Ԍ���,�[,5�Ԍ���,�[,6�ԍ�,�[])
beamin = [0.090,0.090,0.090,0.090,0.080,0.080,0.070,0.070,0.060,0.060,0.035,0.035,0.020,0.025];
%�w��
ply24 = [2 2 2 2 2 2 2 2 1 1 1 1 1 1];
ply40 = [5 5 5 5 5 5 7 7 6 6 6 6 6 6];
%���O�a(m)
beamout = beamin + 0.125*10^(-3)*ply24 + 0.111*10^(-3)*ply40;

%�Ǐ����d��(kg/mm*m/s^2)([1�Ԍ�,2�Ԍ�,3�Ԍ�,4�Ԍ�,5�Ԍ�,6�Ԍ�])
m0beam = [0.35*10^(-3)*g,0.35*10^(-3)*g,0.3*10^(-3)*g,0.25*10^(-3)*g,0.235*10^(-3)*g,0.19*10^(-3)*g,0.08*10^(-3)*g];
%{
���̐�����@�ɂ���ċǏ����d�ʂ͕ς��B19���18��̗����Q�l�ɋǏ����d�ʂ����߂��B
18��̗����@�͗��ǃm�[�g���Q�l�ɂ��Ă��������B
%}
m0wing = [0.75*g 0.7*g 0.7*g 0.65*g 0.5*g 0.4*g,0.3*g];%�Ǐ����d��(kg/m*m/s^2)

%�������O��
E0 = [104.803*10^9 97.009*10^9 97.009*10^9 117.392*10^9 123.059*10^9 123.059*10^9 123.059*10^9 ];


%{
���f�ɂ�����͂��猅�̂���݂��v�Z���Ă���A�z�����ɋy�ڂ��e�����v�Z����B
���̌�œK������B
%}

n = 0;
while n <= N
    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);
    %���C�m���Y��
    [Re] = reynolds(U,chord,nu);
    %�g�͌X��
    [original15m,original15mm] = clslope19(Re,codemax,U,nu);
    %�[���g�͊p
    [alpha0015,alpha0015m] = zerolift19(Re,codemax,U,nu);

    %������͗p�}�p(�x)
    alphas = 0;
    alpha = zeros(1,npartition);
    for i = 1:npartition

      if i <= q1
            alpha(1,i) = (alpha1-alpham)/(q1-0)*(i-0) + alpham + alphas;
      elseif i <= q2
            alpha(1,i) = (alpha2-alpha1)/(q2-q1)*(i-q1) + alpha1 + alphas;
      elseif i <= q3
            alpha(1,i) = (alpha3-alpha2)/(q3-q2)*(i-q2) + alpha2 + alphas;
      elseif i <= q4
            alpha(1,i) = (alpha4-alpha3)/(q4-q3)*(i-q3) + alpha3 + alphas;
      elseif i <= q5
            alpha(1,i) = (alpha5-alpha4)/(q5-q4)*(i-q4) + alpha4 + alphas;
      elseif i <= q6
            alpha(1,i) = (alpha6-alpha5)/(q6-q5)*(i-q5) + alpha5 + alphas;
      end
    end  


    % ���̋Ǐ��d��(kg)
    [ty] = mbeam(m0beam,npartition,q1,q2,q3,q4,q5,q6);
    % ���̋Ǐ��d��(kg)
    [mw] = mwing19(m0wing,npartition,q1,q2,q3,q4,q5,q6,nd);
    %���a
    [bd,bd1,br,br1] = beamd19(beamin,beamsp,npartition,q11,q22,q33,q44,q55,q66,q77,nd);
    %�������O��
    [E] = young19(npartition,E0,q11,q22,q33,q44,q55,q66,q77);
    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    F = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi1 = (sum(D) - sum(dD))/sum(F);

    po1 = (sum(D) - sum(dD))*U;

    %----------------------------
    chord0 = chord0 + cr;

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi2 = (sum(D) - sum(dD))/sum(F);

    po2 = (sum(D) - sum(dD))*U; 

    chord0 = chord0 - cr;

    %--------------------------------------

    chord0 = chord0 - cr;

    if chord0 < chord1
        break
    end    

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi3 = (sum(D) - sum(dD))/sum(F);

    po3 = (sum(D) - sum(dD))*U; 


    chord0 = chord0 + cr;

    %------------------------------------------

    chord1 = chord1 + cr;

    if chord0 < chord1
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi4 = (sum(D) - sum(dD))/sum(F);

    po4 = (sum(D) - sum(dD))*U; 

    chord1 = chord1 - cr;

    %------------------------------------

    chord1 = chord1 - cr;

    if chord1 < chord2
        break
    end    

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi5 = (sum(D) - sum(dD))/sum(F);

    po5 = (sum(D) - sum(dD))*U; 

    chord1 = chord1 + cr;

    %----------------------------------------------------

    chord2 = chord2 + cr;

    if chord1 < chord2
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi6 = (sum(D) - sum(dD))/sum(F);

    po6 = (sum(D) - sum(dD))*U; 

    chord2 = chord2 - cr;

    %-----------------------------------------

    chord2 = chord2 - cr;

    if chord2 < chord3
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi7 = (sum(D) - sum(dD))/sum(F);

    po7 = (sum(D) - sum(dD))*U; 

    chord2 = chord2 + cr;

    %-----------------------------------------

    chord3 = chord3 + cr;

    if chord2 < chord3
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi8 = (sum(D) - sum(dD))/sum(F);

    po8 = (sum(D) - sum(dD))*U; 

    chord3 = chord3 - cr;

    %------------------------------------

    chord3 = chord3 - cr;

    if chord3 < chord4
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi9 = (sum(D) - sum(dD))/sum(F);

    po9 = (sum(D) - sum(dD))*U; 

    chord3 = chord3 + cr;

    %--------------------------------------

    chord4 = chord4 + cr;

    if chord3 < chord4
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi10 = (sum(D) - sum(dD))/sum(F);

    po10 = (sum(D) - sum(dD))*U; 

    chord4 = chord4 - cr;

    %--------------------------------

    chord4 = chord4 - cr;

    if chord4 < chord5
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi11 = (sum(D) - sum(dD))/sum(F);

    po11 = (sum(D) - sum(dD))*U; 

    chord4 = chord4 + cr;

    %-----------------------------------

    chord5 = chord5 + cr;

    if chord4 < chord5
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi12 = (sum(D) - sum(dD))/sum(F);

    po12 = (sum(D) - sum(dD))*U; 

    chord5 = chord5 - cr;

    %-------------------------------------

    chord5 = chord5 - cr;

    if chord5 < chord6
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi13 = (sum(D) - sum(dD))/sum(F);

    po13 = (sum(D) - sum(dD))*U; 

    chord5 = chord5 + cr;

    %--------------------------------------------

    chord6 = chord6 + cr;

    if chord5 < chord6
        break
    end

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi14 = (sum(D) - sum(dD))/sum(F);

    po14 = (sum(D) - sum(dD))*U; 

    chord6 = chord6 - cr;

    %-----------------------------------------------

    chord6 = chord6 - cr;

    %�����_�̍��W(�������̗��}�E���g��x=0,y=0)
    x = zeros(1,npartition);
    for i = 1:npartition
        if i <= q33 
            x(1,i) = nd*i;
        else
            x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
        end
    end

    y = zeros(1,npartition);
    for i = 1:npartition
        if i > sp(1,3)/0.025
            y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
        end
    end

    % chord=������(m)
    [chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

    %�����g��
    [F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
    %�R��
    [drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

    %����݌v�Z

    [Iz] = areainertia19(bd,bd1);
    [M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
    Fm = F.*cos(theta);


    for i = 1:npartition
         if i <= q33

              x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
              y(1,i) = ty(1,i);

         else 
              x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
              y(1,i) = ty(1,i);  

         end
    end

    for i = q33+1:npartition

        x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
        y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

    end

    %�������̏z��
    gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
    %�z���z�v�Z
    [gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
    %���������v�Z
    [X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
    %�z���z�Čv�Z
    [gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
    %�Ǐ��g�͍Čv�Z
    F = rho*U*gammas.*cos(theta)*nd;
    %�Ǐ��U���R��
    dD = rho*downwash.*gammas;

    effi15 = (sum(D) - sum(dD))/sum(F);

    po15 = (sum(D) - sum(dD))*U; 

    chord6 = chord6 + cr;

    %--------------------------------------------

    Effi = [po1 po2 po3 po4 po5 po6 po7 po8 po9 po10 po11 po12 po13 po14 po15];
    MIN = min(Effi);

    if MIN == po1
        break
    end    

    if MIN == po2
        chord0 = chord0 + cr;
    end

    if MIN == po3
        chord0 = chord0 - cr;
    end

    if MIN == po4
        chord1 = chord1 + cr;
    end

    if MIN == po5
        chord1 = chord1 - cr;
    end

    if MIN == po6
        chord2 = chord2 + cr;
    end

    if MIN == po7
        chord2 = chord2 - cr;
    end

    if MIN == po8
        chord3 = chord3 + cr;
    end

    if MIN == po9
        chord3 = chord3 - cr;
    end

    if MIN == po10
        chord4 = chord4 + cr;
    end

    if MIN == po11
        chord4 = chord4 - cr;
    end

    if MIN == po12
        chord5 = chord5 + cr;
    end

    if MIN == po13
        chord5 = chord5 - cr;
    end

    if MIN == po14
        chord6 = chord6 + cr;
    end

    if MIN == po15
        chord6 = chord6 - cr;
    end

    n = n + 1;

    chordtc = [chord0 chord1 chord2 chord3 chord4 chord5 chord6];
    Effic = 1 - (sum(D) - sum(dD))/sum(F);
    dihe = atan(y(1,npartition)/x(1,npartition))*180/pi;

    disp(sum(F))
    disp(sum(D))
    disp(sum(dD))
    disp(chordtc)
    disp(Effic)
    disp(dihe)

    figure(1);
    plot(x,chord);

end

%�����_�̍��W(�������̗��}�E���g��x=0,y=0)
x = zeros(1,npartition);
for i = 1:npartition
    if i <= q33 
        x(1,i) = nd*i;
    else
        x(1,i) = nd*q33 + nd*(i - q33)*cos(thetaC);
    end
end

y = zeros(1,npartition);
for i = 1:npartition
    if i > sp(1,3)/0.025
        y(1,i) = nd*(i - sp(1,3)/0.025)*sin(thetaC);
    end
end

% chord=������(m)
[chord]=CHORD19(chord0,chord1,chord2,chord3,chord4,chord5,chord6,q1,q2,q3,q4,q5,q6,npartition);

%�����g��
[F] = Liftfirst19(original15m,chord,alpha,nd,rho,U,npartition,alpha0015);
%�R��
[drag15,drag15m,D]=drag(Re,codemax,U,nu,alpha,alpham,npartition,chord,rho,nd);

%����݌v�Z

[Iz] = areainertia19(bd,bd1);
[M,theta,ty] = momente19(F,npartition,nd,WW,mw,Iz,E,q33,thetaC,y);
Fm = F.*cos(theta);


for i = 1:npartition
     if i <= q33
         
          x(1,i) = nd*i - tan(theta(1,i))*ty(1,i);
          y(1,i) = ty(1,i);
      
     else 
          x(1,i) = nd*i - tan(theta(1,i)-thetaC)*ty(1,i);
          y(1,i) = ty(1,i);  
     
     end
end

for i = q33+1:npartition
    
    x(1,i) = x(1,q33) + cos(thetaC)*(x(1,i)-x(1,q33)) - sin(thetaC)*(y(1,i)-y(1,q33));
    y(1,i) = y(1,q33) + sin(thetaC)*(x(1,i)-x(1,q33)) + cos(thetaC)*(y(1,i)-y(1,q33));

end

%�������̏z��
gammam = 1/2*original15mm*(alpham-alpha0015m)*U*codemax;
%�z���z�v�Z
[gamma] = gamma19(npartition,original15m,chord,alpha,alpha0015,U);
%���������v�Z
[X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd);
%�z���z�Čv�Z
[gammas] = gammas19(npartition,original15m,chord,alpha,alpha0015,U,downwash);
%�Ǐ��g�͍Čv�Z
F = rho*U*gammas.*cos(theta)*nd;
%�Ǐ��U���R��
dD = rho*downwash.*gammas;


chordtc = [chord0 chord1 chord2 chord3 chord4 chord5 chord6];
Effic = 1 - (sum(D) - sum(dD))/sum(F);
dihe = atan(y(1,npartition)/x(1,npartition))*180/pi;

mac = 2/26.96*0.025*sum(chord.^2);
disp(mac)

disp(sum(F))
disp(sum(D))
disp(sum(dD))
disp(chordtc)
disp(Effic)
disp(dihe)
disp(sum(F)*2/(W*g))
disp((sum(D)-sum(dD))*U*2)
disp(M(1,1))
disp(M(1,WW(1,2)))
disp(ty(1,npartition))

figure(1);
plot(x,chord);
figure(9);
plot(x,downwash);
%theta;
%size(F);
%size(chord);

Cle = 2*rho*U^2*nd*F.*(1./chord);

figure(10);
plot(x,Cle);

ty = ty';
save ty.txt  ty -ascii

