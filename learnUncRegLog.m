function [param_mle,L_mle,grad_L,hessian_L] = learnUncRegLog(X,y)
[N,M]=size(X);
init_param = zeros(1,M+1);
options = optimset('Algorithm','interior-point','MaxIter',100000,'MaxFunEvals',100000,'Display','off');
[param_mle,J,~,~,~,grad_L,hessian_L] = fmincon(@(t)(costFctLogReg(t,X,y)),init_param,[],[],[],[],[],[],[], options);
L_mle=exp(-J);
end


