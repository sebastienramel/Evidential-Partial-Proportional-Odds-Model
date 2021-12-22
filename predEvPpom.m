function [mPred] = predEvPpom(ds,p_ppom_mle,J_ppom_mle,xtst,idnoProp,nreqMax,lambda)
% constants
[n,p] = size(xtst);
q = length(idnoProp);
k = (size(p_ppom_mle,2)+1-p+2*q)/(1+q);

mPred=zeros(n,2^k-1);
for h = 1 : n
    xtarg = xtst(h,:);
    Pishat = predPrPpom(p_ppom_mle,xtarg,idnoProp);
    Pshat = cumsum([0,Pishat]);

    plssPs = zeros(k-1,nreqMax+3);
    Pss = zeros(k-1,nreqMax+3);
    for j=1:k-1
        Phat = Pshat(j+1);
        Ps = [-0.5 0 Phat 1 1.5];
        plPs = [0 0 1 0 0];
        nreq = 0;
        while nreq < nreqMax
            nreq = nreq+1;
            dplPs = [diff(plPs)];
            dPs = [diff(Ps)];    
            areas = abs(dplPs).* dPs;

            % Computation of the score for each known evaluation point of the
            % approximation
            degpt = zeros(length(Ps)-2,1); areasum =  zeros(length(Ps)-2,1);
            for i = 2 : length(Ps)-1
                x_l = Ps(i-1) - Ps(i);
                y_l = plPs(i-1) - plPs(i);
                phi_l = rad2deg(atan(y_l/x_l));

                x_r = Ps(i+1) - Ps(i);
                y_r = plPs(i+1) - plPs(i);
                phi_r = rad2deg(atan(y_r/x_r));

                degpt1 = (180 + phi_l) - phi_r;
                degpt2 = (180 - phi_l) + phi_r;

                degpt(i-1) = min(abs(degpt1),abs(degpt2));
                areasum(i-1) = areas(i-1)+ areas(i);
            end

            % Computation of the value "theta" for the next optimisation round
            [~, imaxdegpl2] = max((180-degpt(2:end-1)).*areasum(2:end-1));
            if areas(imaxdegpl2+1) >  areas(imaxdegpl2+2)
                   P_add =  Ps(imaxdegpl2+1) + 0.5* ( Ps(imaxdegpl2+2) - Ps(imaxdegpl2+1));
            else
                   P_add =  Ps(imaxdegpl2+2) + 0.5* ( Ps(imaxdegpl2+3) - Ps(imaxdegpl2+2));
            end


            % Computation of the contour at the computed value of "theta"
            plLR = cumContPpom(ds,J_ppom_mle,xtarg,idnoProp,P_add,j,lambda);

            % Update and sort known evaluation points
            pl_foo = [plPs,plLR];
            P_foo = [Ps, P_add];
            [Psvals, Psid] = sort(P_foo);
            Ps = Psvals;
            plPs = pl_foo(Psid);

        end
        plssPs(j,:)  = plPs(2:end-1);
        Pss(j,:)  = Ps(2:end-1);

    end
    mPred(h,:) = computePredMass(Pshat,Pss,plssPs);
  


end   

end



