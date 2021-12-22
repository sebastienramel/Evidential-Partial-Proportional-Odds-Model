function [J, grad] = costFctLogReg(param, X, y)

yplatt = zeros(length(y),1);

htheta = 1./ (1+exp(-(param(1)+X*param(2:end)')));


%platt binary lbl transformation
tplus = ( sum(y==1) + 1 )  / ( sum(y==1) + 2); 
tmoins = 1 / (sum(y==0)+2);
yplatt(y==1) = tplus;
yplatt(y==0) = tmoins;

% 
% htheta(htheta==1)= 1-eps;
% htheta(htheta==0)=eps;

%yplatt = y;

% Cost Function : -log L 
J = sum(-yplatt .* log(htheta) - (1 - yplatt) .* log(1 - htheta)) ;
grad = 0;


end
