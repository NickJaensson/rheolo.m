% function to calculate von Mises stresses
function [taud] = von_Mises(svec)

    % unpack the vector
    sxx = svec(1); sxy = svec(2); sxz = svec(3); 
    syy = svec(4); syz = svec(5); szz = svec(6);

    % put the 6 values in a 3x3 symmetrix matrix
    ss = [sxx,sxy,sxz;sxy,syy,syz;sxz,syz,szz];

    % determine von Mises stresses
    sshat = ss - (1/3)*trace(ss)*eye(3);
    taud = sqrt(0.5*(trace(sshat*sshat)));

end