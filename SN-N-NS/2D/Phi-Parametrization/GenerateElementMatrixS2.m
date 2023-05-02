%Compute the matrices obtained in the weak formulation of the 2D Usadel 
%equation in the phi-parametrization in the S2 layer
%S2 corresponds to ksi^2 nabla*(G_S^2 nabla Phi_S)
%M2 corresponds to omega_n*Phi_S*G_S
function [S2elem,M2elem] = GenerateElementMatrixS2(GI,points,weights,omega,ksi,G,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S2elem = ksi^2*abs(GI.Deltas2)/2.*(G.^2).'.*(GI.bs2(:,ind1).*GI.bs2(:,ind2)+GI.cs2(:,ind1).*GI.cs2(:,ind2))*weights.';
    M2elem = abs(GI.Deltas2)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
end