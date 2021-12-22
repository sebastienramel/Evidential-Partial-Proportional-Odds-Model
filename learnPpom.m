function [param_mle,J_mle] = learnPpom(X,idnoProp,y,lambda)
% extraction of constants
q = length(idnoProp);
k = length(unique(y));
p = size(X,2);

% init parameter
alpha_k = [1:k-1].*1e-1;
sigma_vec = ones(1,p-q).*1e-1;
gamma_vec_k = mapparam(repmat([1:(k-1)]'.*1e-1,1,q));
init_param = [alpha_k,sigma_vec,gamma_vec_k];

% linear constrains
Aalph = zeros(k-2,(k-1)+p-q+(k-1)*q); balph = zeros(k-2,1);
%Agamm = zeros((k-2)*q,(k-1)+p-q+(k-1)*q); bgamm = zeros((k-2)*q,1);
for j = 1 : k-2
    Aalph(j,:) = [zeros(1,(j-1)),1,-1,zeros(1,(k-1)-((j-1)+2)),zeros(1,p-q),zeros(1,(k-1)*q)];
    balph(j) = 0;
%     for l = 1 : q
%         Ag = [zeros((j-1),q);zeros(1,(l-1)),1,zeros(1,(q-l));zeros(1,(l-1)),-1,zeros(1,(q-l));zeros((k-1-((j-1)+2)),q)];
%         Agamm((j-1)*q+l,:)=[zeros(1,k-1),zeros(1,p-q),mapparam(Ag)];
%         bgamm((j-1)*q+l) =  0;
%     end
end


% optimisation
options = optimset('Algorithm','interior-point','GradObj', 'on','MaxIter',100000,'MaxFunEvals',100000,'Display','iter');
[param_mle, J_mle] = fmincon(@(param)(costFctPpom(param,X,idnoProp,y,lambda)),init_param,[Aalph],[balph],[],[],[],[],[],options);

end