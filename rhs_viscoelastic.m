function [rhs] = rhs_viscoelastic(cvec,user)

    % define the unit tensor for convenience
    I = eye(3);

    % calculate the velocity gradient tensor
    L = zeros(3);
    if user.flowtype == 1
        L(1,2) = user.rate;
    elseif user.flowtype == 2
        L(1,1) = user.rate; L(2,2) = user.rate;
    elseif user.flowtype == 3
        L(1,1) = user.rate; L(2,2) = -user.rate/2; L(3,3) = -user.rate/2;
    end

    % unpack the vector
    cxx = cvec(1); cxy = cvec(2); cxz = cvec(3); 
    cyy = cvec(4); cyz = cvec(5); czz = cvec(6);

    % put the 6 values in a 3x3 symmetrix matrix
    cc = [cxx,cxy,cxz;cxy,cyy,cyz;cxz,cyz,czz];

    % calculate the upper convected part
    UCDpart = L*cc+cc*transpose(L);

    % calculate the relaxtion part
    if user.model == 1 % UCM
        rlxpart = -1/user.lam*(cc-I);
    elseif user.model == 2 % Giesekus
        rlxpart = -1/user.lam*(cc-I+user.alpha*(cc-I)*(cc-I));
    elseif user.model == 3 % PTTlin
        rlxpart = (-1/user.lam)*(1+user.eps*(trace(cc)-3.0))*(cc-I);
    elseif user.model == 4 % PTTexp
        rlxpart = (-1/user.lam)*(exp(user.eps*(trace(cc)-3.0)))*(cc-I);
    end

    % adapt the relaxation part for the alam models
    if any(user.alam==[2,3])
        stress = stress_viscoelastic_3D(cvec,user);
        taud = von_Mises(stress);
    end

    if user.alam == 1 % elastic
        rlxpart = 0.0;
    elseif user.alam == 2 % Saramito1
        rlxpart = rlxpart*max([0,(taud-user.tauy)/taud]);
    elseif user.alam == 3 % Saramito2
        if taud <= user.tauy
            fac = 0;
        else
            fac = 1/((taud/user.G)*((taud-user.tauy)/user.Kfac)^(-1/user.nexp));
        end
        rlxpart = fac*rlxpart*user.lam;
    end

    % determine the complete solution
    sol = UCDpart+rlxpart;

    % return the 3x3 in a 6 vector
    rhs = [sol(1,1),sol(1,2),sol(1,3),sol(2,2),sol(2,3),sol(3,3)];

end
