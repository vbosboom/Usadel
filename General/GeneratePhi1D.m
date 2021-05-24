%Generates first order lagrangian basis function in a 1D element with 
%endpoints [a,b], length h at point x
function sol = GeneratePhi1D(a,b,h,x)
    sol(1,:) = (b-x)./h;
    sol(2,:) = (x-a)./h;
end