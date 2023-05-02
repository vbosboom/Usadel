%Compute the matrices and Jacobians obtained in the weak formulation of the 2D Usadel 
%equation in the theta,chi-parametrization in the N layer
%S1n corresponds to -(nabla^2 theta_N)
%S2n corresponds to -nabla(sin(theta_N)^2*nabla chi_N)
function [S1nelem, S2nelem, JS23Telem,Jf13Telem,Jf13Celem,Jf23Telem] = GenerateElementMatrixNN(GI,points,weights,E,thetan,Chigrad1,Chigrad2,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1).';
    Phi2 = GeneratePhi(points,ind2).';

    S1nelem = abs(GI.Deltan)/2.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2));
    S2nelem = abs(GI.Deltan)/2.*(sin(thetan).^2).'.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2))*weights.';
    
    JS23Telem = abs(GI.Deltan)/2.*((2*sin(thetan).*cos(thetan)).'.*((GI.bn(:,ind1).*Chigrad1+GI.cn(:,ind1).*Chigrad2).*Phi2))*weights.';
    Jf13Telem = -abs(GI.Deltan)/2.*(cos(thetan).^2-sin(thetan).^2).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1.*Phi2*weights.';
    Jf13Celem = -abs(GI.Deltan)/2.*((2*sin(thetan).*cos(thetan)).'.*(GI.bn(:,ind2).*Chigrad1+GI.cn(:,ind2).*Chigrad2).*Phi1)*weights.';
    Jf23Telem = abs(GI.Deltan)/2.*(1i*E*cos(thetan).').*Phi1.*Phi2*weights.';
end