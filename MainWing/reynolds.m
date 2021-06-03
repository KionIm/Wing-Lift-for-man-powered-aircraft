function [Re] = reynolds(U,chord,nu)

%U:機速、cy:コード長、nu:動粘性係数
    
    Re = U*chord/nu;
    
end