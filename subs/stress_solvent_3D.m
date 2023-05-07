% calculate the solvent stress
function [tauvec] = stress_solvent_3D(vemodel,rate,flowtype)

    % calculate the velocity gradient tensor
    L = fill_L(vemodel,rate,flowtype);

    sol = vemodel.eta_s * (L+transpose(L));

    % return the 3x3 in a 6 vector
    tauvec = [sol(1,1),sol(1,2),sol(1,3),sol(2,2),sol(2,3),sol(3,3)];
end