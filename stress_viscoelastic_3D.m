% calculate the viscoelastic stress
function [tauvec] = stress_viscoelastic_3D(cvec)

    global G

    tauvec(1) = G*(cvec(1)-1);
    tauvec(2) = G*(cvec(2));
    tauvec(3) = G*(cvec(3));
    tauvec(4) = G*(cvec(4)-1);
    tauvec(5) = G*(cvec(5));
    tauvec(6) = G*(cvec(6)-1);

end