% calculate the viscoelastic stress
function [tauvec] = stress_viscoelastic_3D(cvec,vemodel)

    tauvec(1) = vemodel.G*(cvec(1)-1);
    tauvec(2) = vemodel.G*(cvec(2));
    tauvec(3) = vemodel.G*(cvec(3));
    tauvec(4) = vemodel.G*(cvec(4)-1);
    tauvec(5) = vemodel.G*(cvec(5));
    tauvec(6) = vemodel.G*(cvec(6)-1);

end