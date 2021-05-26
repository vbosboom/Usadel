function [S3melem,S3pelem,M3melem,M3pelem] = GenerateVecElementMatrixF(GI,points,weights,omega_T,G,ind1,ind2)

    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S3melem = abs(GI.Deltan)/2.*(G.^2).'.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2))*weights.';
    S3pelem = abs(GI.Deltan)/2.*(G.^2).'.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2))*weights.';
    M3melem = abs(GI.Deltan)/2.*(omega_T*G).'.*Phi1.*Phi2*weights.';
    M3pelem = abs(GI.Deltan)/2.*(omega_T*G).'.*Phi1.*Phi2*weights.';
end