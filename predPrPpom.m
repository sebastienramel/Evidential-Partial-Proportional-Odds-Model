function [pi,deci] = predPrPpom(p_ppom_mle,Xtst,idnoProp)
[n, p] = size(Xtst);
q = length(idnoProp);
Xnpo = Xtst(:,ismember(1:p,idnoProp)); % features satisfying proportionality assumption
Xpo = Xtst(:,~ismember(1:p,idnoProp)); % features violating proportionality assumption
k = (size(p_ppom_mle,2)+1-p+2*q)/(1+q); % # of classe labels

% map parameter
alpha_k = p_ppom_mle(1:k-1);
sigma_vec = p_ppom_mle(k:k+p-q-1);
gamma_vec_k = demapparam(p_ppom_mle(k+p-q:end),q);


% Computation 
Sigma = [alpha_k',-repmat(sigma_vec,k-1,1),gamma_vec_k];
Z = (Sigma*[ones(n,1),Xpo,Xnpo]')';
h = 1./(1+exp(-Z));
hh = [zeros(n,1),h,ones(n,1)];%[h(z0),h(z1),...,h(z(k-1)),h(zk)]
pi = diff(hh')';
[~,deci] = max(pi,[],2);

end

