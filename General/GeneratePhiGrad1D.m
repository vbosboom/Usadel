%Evalutes the gradient of the two linear basis function on a
%1D element with length h.
function sol = GeneratePhiGrad1D(h)
    sol(1,1) = -1/h;
    sol(2,1) = 1/h;
end