%Evaluates the two linear basis function on a 1D element with 
%endpoints [a,b] and length h at the point x
function sol = GeneratePhi1D(a,b,h,x)
    sol(1,:) = (b-x)./h;
    sol(2,:) = (x-a)./h;
end