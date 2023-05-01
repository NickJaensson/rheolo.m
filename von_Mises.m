% function to calculate von Mises stresses
function [taud] = von_Mises(stress)

    stresshat = stress - (1/3)*trace(stress)*eye(3);
    taud = sqrt(0.5*(trace(stresshat*stresshat)));

end