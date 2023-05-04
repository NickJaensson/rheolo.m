% calculate the solvent stress
function [tauvec] = stress_solvent_3D(vemodel)

    % calculate the velocity gradient tensor
    L = zeros(3);
    if vemodel.flowtype == 1
        L(1,2) = vemodel.rate;
    elseif vemodel.flowtype == 2
        L(1,1) = vemodel.rate; L(2,2) = vemodel.rate;
    elseif vemodel.flowtype == 3
        L(1,1) = vemodel.rate; L(2,2) = -vemodel.rate/2; L(3,3) = -vemodel.rate/2;
    end

    sol = vemodel.eta_s * (L+transpose(L));

    % return the 3x3 in a 6 vector
    tauvec = [sol(1,1),sol(1,2),sol(1,3),sol(2,2),sol(2,3),sol(3,3)];
end