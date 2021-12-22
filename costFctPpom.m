function [J,grad] = costFctPpom(p_ppom, X, idnoProp, y, lambda)
% constants
[n, p] = size(X);
q = length(idnoProp);
k = (size(p_ppom,2)+1-p+2*q)/(1+q); % # of classe labels

% dummy variable label
dummylbls = eye(k);
t = dummylbls(y,:);

% /!\ regularization with platt labels transformation (lambda=0) /!\
lambda=0;
tt = zeros(size(t));teye = eye(k);
for kk = 1:k
    tplus = ( sum(y==kk) + 1 )  / ( sum(y==kk) + k); 
    tmoins = 1 / (sum(y==kk)+k);
    aa = teye(kk,:);    bb = teye(kk,:);
    bb(aa==1) = tplus; bb(aa==0)=tmoins;
    tt(bi2de(t)==bi2de(teye(kk,:)),:) = repmat(bb,sum(bi2de(t)==bi2de(teye(kk,:))),1);
end
t = tt;


% map features
Xnpo = X(:,ismember(1:p,idnoProp)); % features satisfying proportionality assumption
Xpo = X(:,~ismember(1:p,idnoProp)); % features violating proportionality assumption

% map parameters
alpha_k = p_ppom(1:k-1);
sigma_vec = p_ppom(k:k+p-q-1);
gamma_vec_k = demapparam(p_ppom(k+p-q:end),q);

% Computation 
Sigma = [alpha_k',-repmat(sigma_vec,k-1,1),gamma_vec_k];
Z = (Sigma*[ones(n,1),Xpo,Xnpo]')';
h = 1./(1+exp(-Z));
hh = [zeros(n,1),h,ones(n,1)];%[h(z0),h(z1),...,h(z(k-1)),h(zk)]
pi = diff(hh')';
dhdze = [zeros(n,1),(1./(1+exp(-Z))).*(1./(1+exp(Z))),zeros(n,1)];%[h(z0)h(-z0),h(z1)h(-z1),...,h(zk)h(-zk)]
diffdhdze = diff(dhdze')';

temp_p_ppom = [zeros(1,k-1),p_ppom(k:end)];

% Cost Function : -log L 
J = -sum(sum(t .* log(pi),2)) +  lambda/2*sum(temp_p_ppom.^2) ;

% Gradient 
grad = zeros(size(p_ppom));
for j=1:k-1
    grad(j)= -sum(dhdze(:,j+1).*(t(:,j)./pi(:,j) - t(:,j+1)./pi(:,j+1))) + lambda*temp_p_ppom(j);  
end
for l=1:p-q
    grad(k-1+l)= sum(sum(diffdhdze .* t./pi .* repmat(Xpo(:,l),1,k), 2)) + lambda*temp_p_ppom(j+l);
end
for jj=1:k-1
    for ll=1:q
        grad((k-1)+(p-q)+(jj-1)*q+ll) = -sum(dhdze(:,jj+1).* (Xnpo(:,ll).*(t(:,jj)./pi(:,jj) - t(:,jj+1)./pi(:,jj+1)))) + lambda*temp_p_ppom((k-1)+(p-q)+(jj-1)*q+ll);
    end
end
a= 0;