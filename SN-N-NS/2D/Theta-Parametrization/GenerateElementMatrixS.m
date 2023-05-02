%Compute the matrices obtained in the weak formulation of the 2D Usadel 
%equation in the theta,chi-parametrization in the S layer
%S1 corresponds to -ksi^2*(nabla^2 theta_S)
%S2 corresponds to -ksi^2*nabla(sin(theta_S)^2*nabla chi_S)
function [S1elem,S2elem] = GenerateElementMatrixS(Delta,b,c,weights,ksi,Theta,ind1,ind2)
    S1elem = ksi^2*abs(Delta)/2.*(b(:,ind1).*b(:,ind2)+c(:,ind1).*c(:,ind2));
    S2elem = ksi^2*abs(Delta)/2.*(sin(Theta).^2).'.*(b(:,ind1).*b(:,ind2)+c(:,ind1).*c(:,ind2))*weights.';
end