function [S1elem,S2elem] = GenerateVecElementMatrix(weights,Theta,GradPhi,h,ind1,ind2)
    S1elem = h/2.*(GradPhi(:,ind1).*GradPhi(:,ind2))*2;
    S2elem = h/2.*(sin(Theta).^2.*GradPhi(:,ind1).*GradPhi(:,ind2))*weights;
end