function [rhs] = rhs_viscoelastic(cvec,L,vemodel)

    % define the unit tensor for convenience
    I = eye(3);

    % unpack the vector
    cxx = cvec(1); cxy = cvec(2); cxz = cvec(3); 
    cyy = cvec(4); cyz = cvec(5); czz = cvec(6);

    % put the 6 values in a 3x3 symmetrix matrix
    cc = [cxx,cxy,cxz;cxy,cyy,cyz;cxz,cyz,czz];

    % calculate the upper convected part
    UCDpart = L*cc+cc*transpose(L);

    % calculate the relaxtion part
    if vemodel.model == 1 % UCM
        rlxpart = -1/vemodel.lam*(cc-I);
    elseif vemodel.model == 2 % Giesekus
        rlxpart = -1/vemodel.lam*(cc-I+vemodel.alpha*(cc-I)*(cc-I));
    elseif vemodel.model == 3 % PTTlin
        rlxpart = (-1/vemodel.lam)*(1+vemodel.eps*(trace(cc)-3.0))*(cc-I);
    elseif vemodel.model == 4 % PTTexp
        rlxpart = (-1/vemodel.lam)*(exp(vemodel.eps*(trace(cc)-3.0)))*(cc-I);
    end

    % adapt the relaxation part for the alam models
    if any(vemodel.alam==[2,3])
        stress = stress_viscoelastic_3D(cvec,vemodel);
        taud = von_Mises(stress);
    end

    if vemodel.alam == 1 % elastic
        rlxpart = 0.0;
    elseif vemodel.alam == 2 % Saramito1
        rlxpart = rlxpart*max([0,(taud-vemodel.tauy)/taud]);
    elseif vemodel.alam == 3 % Saramito2
        if taud <= vemodel.tauy
            fac = 0;
        else
            fac = 1/((taud/vemodel.G)*((taud-vemodel.tauy)/vemodel.Kfac)^(-1/vemodel.nexp));
        end
        rlxpart = fac*rlxpart*vemodel.lam;
    end

    % determine the complete solution
    sol = UCDpart+rlxpart;

    % return the 3x3 in a 6 vector
    rhs = [sol(1,1),sol(1,2),sol(1,3),sol(2,2),sol(2,3),sol(3,3)];

end
