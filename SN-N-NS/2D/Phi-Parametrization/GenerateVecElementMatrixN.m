function [S3elem,M3elem] = GenerateVecElementMatrixN(GI,points,weights,omega,ksi,G,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S3elem = abs(GI.Deltan)/2.*(G.^2).'.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2))*weights.';
    M3elem = abs(GI.Deltan)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
end