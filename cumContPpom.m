function [plPj] = cumContPpom(ds,J_ppom_mle,xtarg,idnoProp,Pjj,j,lambda)

% extraction of constants
p = size(xtarg,2);
q = length(idnoProp);
k = length(unique(ds(:,end)));

% map features
Xtnpo = xtarg(:,ismember(1:p,idnoProp)); % features satisfying proportionality assumption
Xtpo = xtarg(:,~ismember(1:p,idnoProp)); % features violating proportionality assumption

% init parameter
sigma_vec = ones(1,p-q).*1e-1;
gamma_vec_k = mapparam(repmat([1:(k-1)]'.*1e-1,1,q));
gamm = repmat([1:(k-1)]'.*1e-1,1,q);

if q>0
    alphaj0(j) = log(1/Pjj-1)+ [sigma_vec,-gamm(j,:)]*[Xtpo,Xtnpo]';            
else
    alphaj0(j) = log(1/Pjj-1)+ [sigma_vec]*[Xtpo]';            
end
for a = (j-1) : - 1 : 1
    alphaj0(a) = alphaj0(a+1)-1;
end
for b = (j+1) : 1 : k-1
    alphaj0(b) = alphaj0(b-1)+1;
end
p_ppom_0 = [alphaj0,sigma_vec,gamma_vec_k];

%linear constrains for the class depedent parameters
Aalph = zeros(k-2,(k-1)+p-q+(k-1)*q); balph = zeros(k-2,1);
%Agamm = zeros((k-2)*q,(k-1)+p-q+(k-1)*q); bgamm = zeros((k-2)*q,1);
for r = 1 : k-2
    Aalph(r,:) = [zeros(1,(r-1)),1,-1,zeros(1,(k-1)-((r-1)+2)),zeros(1,p-q),zeros(1,(k-1)*q)];
    balph(r) = 0;
%     for l = 1 : q
%         Ag = [zeros((r-1),q);zeros(1,(l-1)),1,zeros(1,(q-l));zeros(1,(l-1)),-1,zeros(1,(q-l));zeros((k-1-((r-1)+2)),q)];
%         Agamm((r-1)*q+l,:)=[zeros(1,k-1),zeros(1,p-q),mapparam(Ag)];
%         bgamm((r-1)*q+l) =  0;
%     end
end

% equality constrain for evidential forecasting
subAAeq = [zeros(j-1,q);ones(1,q).*Xtnpo;zeros(k-j-1,q)];
AAeq=[zeros(1,j-1),1,zeros(1,k-j-1),-Xtpo,mapparam(subAAeq)];
bbeq= -log(1/Pjj-1);

options = optimset('Algorithm','interior-point','GradObj', 'on','MaxIter',100000,'MaxFunEvals',100000,'Display','off');
[betaoptim ,J_optim]= fmincon(@(beta)(costFctPpom(beta, ds(:,1:end-1), idnoProp, ds(:,end), lambda)),p_ppom_0,[Aalph],[balph],AAeq,bbeq,[],[],[],options);
plPj = exp(J_ppom_mle-J_optim);

end

