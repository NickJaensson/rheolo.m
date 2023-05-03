% calculate the viscoelastic stress
function [tauvec] = stress_viscoelastic_3D(cvec,user)

    tauvec(1) = user.G*(cvec(1)-1);
    tauvec(2) = user.G*(cvec(2));
    tauvec(3) = user.G*(cvec(3));
    tauvec(4) = user.G*(cvec(4)-1);
    tauvec(5) = user.G*(cvec(5));
    tauvec(6) = user.G*(cvec(6)-1);

end