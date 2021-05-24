%Generates gradient of a first order Lagrangian basis function on an
%element with length h.
function sol = GeneratePhiGrad1D(h)
    sol(1,1) = -1/h;
    sol(2,1) = 1/h;
end