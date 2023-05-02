%Compute the matrices obtained in the weak formulation of the 2D Usadel 
%equation in the theta,chi-parametrization in the N layer
%S1n corresponds to -*(nabla^2 theta_N)
%S2n corresponds to -*nabla(sin(theta_N)^2*nabla chi_N)
function [S1nelem, S2nelem] = GenerateElementMatrixN(GI,weights,thetan,ind1,ind2)
    S1nelem = abs(GI.Deltan)/2.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2));
    S2nelem = abs(GI.Deltan)/2.*(sin(thetan).^2).'.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2))*weights.';
end