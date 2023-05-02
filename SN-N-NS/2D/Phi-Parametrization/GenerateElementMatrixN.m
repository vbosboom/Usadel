%Compute the matrices obtained in the weak formulation of the 2D Usadel 
%equation in the phi-parametrization in the N layer
%S3 corresponds to nabla*(G_N^2 nabla Phi_N)
%M3 corresponds to omega_n*Phi_N*G_N
function [S3elem,M3elem] = GenerateElementMatrixN(GI,points,weights,omega,G,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S3elem = abs(GI.Deltan)/2.*(G.^2).'.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2))*weights.';
    M3elem = abs(GI.Deltan)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
end