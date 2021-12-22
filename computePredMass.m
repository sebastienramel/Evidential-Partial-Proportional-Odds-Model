function [mPred] = computePredMass(Pshat,Pss,plssPs)
K = length(Pshat)-1;
delt = diff(Pshat); epsi=1e-4;
idpinotnull = find(delt>epsi); 
idpinull = find(delt<epsi);
k=length(idpinotnull);
idcomputeCont = idpinotnull;


if length(idcomputeCont)==1
    mPred = zeros(1,2^K-1);
    if idcomputeCont==K
        mPred(bi2de([zeros(1,K-1),1]))=1;
    elseif idcomputeCont==1
        mPred(bi2de([1,zeros(1,K-1)]))=1;
    end
else
    idcomputeCont(idcomputeCont==K)=[];
    TabSig = zeros(k,k);
    lowEnvSig = zeros(k,k); 
    uppEnvSig = zeros(k,k); 
    lowEnvSig(1,:)=idpinotnull;
    uppEnvSig(1,:)=idpinotnull;
    for j = 1 : k
        TabSig(1,j)=bi2de([zeros(1,lowEnvSig(1,j)-1),ones(1,uppEnvSig(1,j)-(lowEnvSig(1,j)-1)),zeros(1,k-(uppEnvSig(1,j)+1))]);
    end
    for l = 2 : k 
        for j = 1 : k-(l-1)
            lowEnvSig(l,j)=lowEnvSig(l-1,j);
            uppEnvSig(l,j)=uppEnvSig(l-1,j+1);
            TabSig(l,j)=bi2de([zeros(1,lowEnvSig(l,j)-1),ones(1,uppEnvSig(l,j)-(lowEnvSig(l,j)-1)),zeros(1,k-(uppEnvSig(l,j)+1))]);
        end
    end


    TabOme = zeros(k,k);
    lowEnvOme = zeros(k,k);
    uppEnvOme = zeros(k,k);
    lowEnvOme(1,:)=1:length(idpinotnull);
    uppEnvOme(1,:)=1:length(idpinotnull);
    for j = 1 : k
        TabOme(1,j)=bi2de([zeros(1,lowEnvOme(1,j)-1),ones(1,uppEnvOme(1,j)-(lowEnvOme(1,j)-1)),zeros(1,k-(uppEnvOme(1,j)+1))]);
    end
    for l = 2 : k 
        for j = 1 : k-(l-1)
            lowEnvOme(l,j)=lowEnvOme(l-1,j);
            uppEnvOme(l,j)=uppEnvOme(l-1,j+1);
            TabOme(l,j)=bi2de([zeros(1,lowEnvOme(l,j)-1),ones(1,uppEnvOme(l,j)-(lowEnvOme(l,j)-1)),zeros(1,k-(uppEnvOme(l,j)+1))]);
        end
    end


    PshatUniq = Pshat;
    if ~isempty(idpinull)
        PshatUniq(idpinull+1)=[];
    end 
    [mPredOme] = computeAreaUniqCont(PshatUniq,Pss(idcomputeCont,:),plssPs(idcomputeCont,:));


    mPredSig = zeros(1,2^K-1);
    for l = 1 : k 
        for j = 1 : k-(l-1)
            mPredSig(1,TabSig(l,j))=mPredOme(1,TabOme(l,j));
        end
    end

    mPred = mPredSig;


end











end

