function [rhs] = rhs_viscoelastic(cvec)

    global model flowtype rate lam alpha eps G alam tauy Kfac nexp

    % convenient definition
    I = eye(3);

    % calculate the velocity gradient tensor
    L = zeros(3);
    if flowtype == 1
        L(1,2) = rate;
    elseif flowtype == 2
        L(1,1) = rate; L(2,2) = rate;
    elseif flowtype == 3
        L(1,1) = rate; L(2,2) = -rate/2; L(3,3) = -rate/2;
    end

    % unpack the vector
    cxx = cvec(1); cxy = cvec(2); cxz = cvec(3); 
    cyy = cvec(4); cyz = cvec(5); czz = cvec(6);

    % put the 6 values in a 3x3 symmetrix matrix
    cc = [cxx,cxy,cxz;cxy,cyy,cyz;cxz,cyz,czz];

    % calculate the upper convected part
    UCDpart = L*cc+cc*transpose(L);

    % calculate the relaxtion part
    if model == 1 % UCM
        rlxpart = -1/lam*(cc-I);
    elseif model == 2 % Giesekus
        rlxpart = -1/lam*(cc-I+alpha*(cc-I)*(cc-I));
    elseif model == 3 % PTTlin
        rlxpart = (-1/lam)*(1+eps*(trace(cc)-3.0))*(cc-I);
    elseif model == 4 % PTTexp
        rlxpart = (-1/lam)*(exp(eps*(trace(cc)-3.0)))*(cc-I);
    end

    % adapt the relaxation part for the alam models
    if any(alam==[1,2,3])
        taud = von_Mises(G*(cc-I));
    end

    if alam == 1 % elastic
        rlxpart = 0.0;
    elseif alam == 2 % Saramito1
        rlxpart = rlxpart*max([0,(taud-tauy)/taud]);
    elseif alam == 3 % Saramito2
        if taud <= tauy
            fac = 0;
        else
            fac = 1/((taud/G)*((taud-tauy)/Kfac)^(-1/nexp));
        end
        rlxpart = fac*rlxpart*lam;
    end

    % determine the complete solution
    sol = UCDpart+rlxpart;

    % return the 3x3 in a 6 vector
    rhs = [sol(1,1),sol(1,2),sol(1,3),sol(2,2),sol(2,3),sol(3,3)];

end
