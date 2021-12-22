rng(1); N=30; % number of instances per class

% Generation of training instances with 3 ordinal classes
mu1 =[0, 0] ; mu2=[3,3]; mu3 = [6,6];
S1=[2 0;0 2]; S2=[2 0;0 2]; S3 = [2 0;0 2];
Xfeat = [mvnrnd(mu1,S1,N);mvnrnd(mu2,S2,N);mvnrnd(mu3,S3,N)];
ds= [Xfeat,[ones(N,1);2*ones(N,1);3*ones(N,1)]];
ds = ds(randperm(3*N),:);

% Learn the partial proportional odds model by finding its MLE
pvalues = brant(ds);
idnoProp = find(pvalues(2:end)<0.05);
[p_ppom_mle, J_ppom_mle]=learnPpom(ds(:,1:end-1),idnoProp,ds(:,end),0);

% Predict the class of a test instance (predictive probability and mass function)
xt = [0,8];
[pit,~] = predPrPpom(p_ppom_mle,xt,idnoProp);
[mt]=predEvPpom(ds,p_ppom_mle,J_ppom_mle,xt,idnoProp,20,0); 

% Decisions based on bayesian and interval dominance strategies with 0/1 cost
[~, decibayes] = max(pit,[],2);
[deciid, ~] = intervalDominance(mt);  