function [mPred] = computeAreaUniqCont(Pshat,Pss,plssPs)
K = length(Pshat)-1;
tempxq=0:1e-5:1;
xq = sort([Pshat(2:end-1),tempxq]);


slj = zeros(K+1,K+1);
plinterpJ = zeros(K-1,length(xq));

% Interpolation of contour functions to find sadle points
%figure;hold on;
idq = find(xq==Pshat(1));
slj(1,1)=idq;
for j=1:K-1
    plinterpJ(j,:)=interp1(Pss(j,:),plssPs(j,:),xq);
    %plot(xq,plinterpJ(j,:),'-k');   
    
    idq = find(xq==Pshat(j+1));
    %plot(xq(idq),1,'xk');
    slj(1,j+1) = idq;
end
idq = find(xq==1);
slj(1,K+1)=idq;

for l = 1:K-2 % level of crossing contour
    slj(l+1,1)=1;
    for j = 1 : K-1-l % id of crossing point
        jj = j+l;     
        [idq,plcross] = findCrossingP(xq,plinterpJ(jj-l,:),plinterpJ(jj,:));
        if idq < slj(l,j+1)
            idq = slj(l,j+1);
        end
        if idq > slj(l,j+2)
            idq = slj(l,j+1);
        end
        slj(l+1,j+1)=idq;
%         switch l
%             case 1
%                 plot(xq(idq),plcross,'ok')
%             case 2  
%                 plot(xq(idq),plcross,'^k')
%             case 3
%                 plot(xq(idq),plcross,'sk')
%             case 4
%                 plot(xq(idq),plcross,'dk')
%             otherwise
%                 plot(xq(idq),plcross,'|k')                
%         end
    end
    idq = find(xq==1);
    slj(l+1,K-l+1)=idq;
end
slj(K,1)=1;slj(K,2)=idq;
slj(K+1,1)=1;


plinterpJext = [[1,zeros(1,size(plinterpJ,2)-1)];plinterpJ;[zeros(1,size(plinterpJ,2)-1),1]];
mPred = zeros(1,2^K-1);
for impl=1:2^K-1
    idj = find(de2bi(impl,K)==1);    
    if sum(de2bi(impl,K))==1 % singletons 
        A_0J =  Pshat(idj+1)-Pshat(idj);
        Lp = [plinterpJext(idj,slj(1,idj):slj(2,idj)),plinterpJext(idj+1,slj(2,idj)+1:slj(1,idj+1))];
        B_0J = trapz(xq(slj(1,idj):slj(1,idj+1)),Lp);
        mPred(1,impl)= A_0J - B_0J;       
    else 
        cnd=(diff(find(de2bi(impl,K)==1)) == 1); % intervals
        if all(cnd)
            j=idj(end);
            nj = length(idj);i = nj-1;
            Up = [plinterpJext(j,slj(nj,j-i):slj(nj-1,j-i+1)), plinterpJext(j-i+1,slj(nj-1,j-i+1)+1:slj(nj,j-i+1))];            
            A_IJ = trapz(xq(slj(nj,j-i):slj(nj,j-i+1)),Up);
            
            Lp= [plinterpJext(j-i,slj(nj,j-i):slj(nj+1,j-i)), plinterpJext(j+1,slj(nj+1,j-i)+1:slj(nj,j-i+1))];
            B_IJ = trapz(xq(slj(nj,j-i):slj(nj,j-i+1)),Lp);
            mPred(1,impl)= A_IJ - B_IJ;
         end
    end
end


probmPred = mPred./sum(mPred);
if sum(mPred)>1
    deltpos = sum(mPred)-1;
    mPred = mPred - probmPred * deltpos;
elseif sum(mPred)<1
    deltneg = 1-sum(mPred);
    mPred = mPred + probmPred * deltneg;
end


end

