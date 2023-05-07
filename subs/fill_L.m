% function to fill velocity gradient tensor
function [L] = fill_L(vemodel,rate,flowtype)

    % calculate the velocity gradient tensor
    L = zeros(3);
    if flowtype == 1
        L(1,2) = rate;
    elseif flowtype == 2
        L(1,1) = rate; L(2,2) = rate;
    elseif flowtype == 3
        L(1,1) = rate; L(2,2) = -rate/2; L(3,3) = -rate/2;
    end

end