function [mappedparameter] = mapparam(demappedparameter)
[k,m] = size(demappedparameter);
mappedparameter = zeros(1,k*m);
for i = 1 : k
    for ii = 1 : m
            mappedparameter(1,(i-1)*m+ii) = demappedparameter(i,ii);
    end
end
%mappedparameter
end

