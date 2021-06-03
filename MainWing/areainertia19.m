function [Iz] = areainertia19(bd,bd1)

%bd:Œ…ŠOŒaAbd1:Œ…“àŒa
%Iz:’f–Ê“ñŸƒ‚[ƒƒ“ƒg

Iz = (pi/64)*(bd.^4-bd1.^4);

end
