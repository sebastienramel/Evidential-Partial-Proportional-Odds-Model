function [demappedparameter] = demapparam(mappedparameter,m)
if ~isempty(mappedparameter)
    kkk = ceil( length(mappedparameter)/m);
    iii =  length(mappedparameter)-(kkk-1)*m;  
    demappedparameter = zeros(kkk,iii);
    for i = 1 : length(mappedparameter)
        k=ceil(i/m);
        ii = i-(k-1)*m;
        demappedparameter(k,ii) = mappedparameter(i);
    end
else
    demappedparameter = [];
end
