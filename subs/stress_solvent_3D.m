% calculate the solvent stress
function [tauvec] = stress_solvent_3D(user)

    % calculate the velocity gradient tensor
    L = zeros(3);
    if user.flowtype == 1
        L(1,2) = user.rate;
    elseif user.flowtype == 2
        L(1,1) = user.rate; L(2,2) = user.rate;
    elseif user.flowtype == 3
        L(1,1) = user.rate; L(2,2) = -user.rate/2; L(3,3) = -user.rate/2;
    end

    sol = user.eta_s * (L+transpose(L));

    % return the 3x3 in a 6 vector
    tauvec = [sol(1,1),sol(1,2),sol(1,3),sol(2,2),sol(2,3),sol(3,3)];
end