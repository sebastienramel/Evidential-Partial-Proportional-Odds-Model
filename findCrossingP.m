function [Idcross, Pcross] = findCrossingP(xq, plPi,plPj)
% with i<j

Pihat = xq(find(plPi==1));
Pjhat = xq(find(plPj==1));

C = find( (plPi-plPj)<0 & xq<Pjhat & xq>Pihat);
if isempty(C)
    Idcross=0;
else
    Idcross = C(1);
end

if Idcross >0
Pcross = plPi(Idcross);
else
Pcross=-1;
end
% a = 0;
end

