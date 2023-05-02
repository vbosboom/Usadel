%Compute the matrices obtained in the weak formulation of the 2D Usadel 
%equation in the phi-parametrization in the S1 material
%S1 corresponds to ksi^2 nabla*(G_S^2 nabla Phi_S)
%M1 corresponds to omega_n*Phi_S*G_S
function [S1elem,M1elem] = GenerateElementMatrixS1(GI,points,weights,omega,ksi,G,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S1elem = ksi^2*abs(GI.Deltas1)/2.*(G.^2).'.*(GI.bs1(:,ind1).*GI.bs1(:,ind2)+GI.cs1(:,ind1).*GI.cs1(:,ind2))*weights.';
    M1elem = abs(GI.Deltas1)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
end