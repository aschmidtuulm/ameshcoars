function [coordinates,elements4,elements3] = ...
    QrefineRG(coordinates,elements3,elements4,marked3,marked4)

[coordinates,elements4,marked,irregular] = ...
    recoarseedges_tri(coordinates,elements3,elements4,marked3,marked4)
[coordinates,elements4,dirichlet,neumann] = ...
    QrefineR(coordinates,elements4,dirichlet,neumann,marked)
[coordinates,elements4,elements3] = ...
    regularizeedges_tri(coordinates,elements4,irregular)
end