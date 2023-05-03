% function to calculate von Mises stresses
function [taud] = von_Mises(svec)

    taud = sqrt ( ( ( svec(1) - svec(4) ) ^ 2 + ...
                    ( svec(4) - svec(6) ) ^ 2 + ...
                    ( svec(6) - svec(1) ) ^ 2 ) / 6 + ...
                      svec(2) ^ 2 + svec(3) ^ 2 + svec(5) ^ 2 );

end