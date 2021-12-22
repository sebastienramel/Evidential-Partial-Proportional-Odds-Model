function [pvalues] = brant(ds)
X = ds(:,1:end-1);
y = ds(:,end);
k = length(unique(y));
[n,p] = size(X);

pj_mle = zeros(k-1,p+1);
Lj_mle = zeros(k-1);
pihat = zeros(n,k-1);
cov_pj_mle = zeros(k-1,p+1,p+1);
for j = 1 : k-1
    z = y>j;
    [pj_mle(j,:),Lj_mle(j),~,H_Lj]=learnUncRegLog(X,z);
    pihat(:,j) = 1./ (1+exp(-(pj_mle(j,1)+X*pj_mle(j,2:end)')));
    cov_pj_mle(j,:,:) = inv(H_Lj);
end

Xplus = [ones(n,1), X];
for m = 1 : (k-2)
    for l = (m+1) : (k-1)
        Wml = diag(pihat(:,l)-pihat(:,m).*pihat(:,l));
        Wm  = diag(pihat(:,m)-pihat(:,m).*pihat(:,m));
        Wl  = diag(pihat(:,l)-pihat(:,l).*pihat(:,l));
        A = inv(Xplus'*Wm*Xplus)*(Xplus'*Wml*Xplus)*inv(Xplus'*Wl*Xplus);
        A(1,:)=[];A(:,1)=[];
        varParam(((m-1)*p+1):(m*p),((l-1)*p+1):(l*p)) = A;
        varParam(((l-1)*p+1):(l*p),((m-1)*p+1):(m*p)) = varParam(((m-1)*p+1):(m*p),((l-1)*p+1):(l*p)) ; 
    end
end

p_concat = [];
for m = 1 : (k-1)
    pmle = pj_mle(m,:); pmle(:,1)=[];
    p_concat = [p_concat; pmle'];
end

for m = 1 : (k-1)
    covparam = squeeze(cov_pj_mle(m,:,:));
    covparam(1,:)=[]; covparam(:,1)=[];
    varParam(((m-1)*p+1):(m*p),((m-1)*p+1):(m*p)) = covparam;
end

I = diag(ones(1,p));
E0 = diag(zeros(1,p));
for i = 1 : k-2
    for j = 1 : k-1
        if j==1
            temp=I;
        elseif j==i+1
            temp=[temp,-I];
        else
            temp=[temp,E0];
        end
    end
    if i==1
        D=temp;
    else
        D=[D;temp];
    end
end
X2 = (D*p_concat)'*inv(D*varParam*D')*(D*p_concat);
df = (k-2)*p;

for kk = 1 : p
    seqs = kk:p:p*(k-1);
    Dseqs = D(:,seqs);
    Dseqs = Dseqs(find(any(Dseqs,2)),:);
    if ~isempty(Dseqs)
        X2 = [X2;(Dseqs*p_concat(seqs))'*inv(Dseqs*varParam(seqs,seqs)*Dseqs')*(Dseqs*p_concat(seqs))];
    else
        X2 = [X2;(Dseqs*p_concat(seqs))'*inv(Dseqs*varParam(seqs,seqs)*(Dseqs')')*(Dseqs*p_concat(seqs))];
    end
    df = [df; (k-2)];
end

pvalues = chi2cdf(X2,df,'upper');
%pvalues
%X2


end

