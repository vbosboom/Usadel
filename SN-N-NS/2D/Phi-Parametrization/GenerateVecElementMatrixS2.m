function [S2elem,M2elem] = GenerateVecElementMatrixS2(GI,points,weights,omega,ksi,G,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S2elem = ksi^2*abs(GI.Deltas2)/2.*(G.^2).'.*(GI.bs2(:,ind1).*GI.bs2(:,ind2)+GI.cs2(:,ind1).*GI.cs2(:,ind2))*weights.';
    M2elem = abs(GI.Deltas2)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
end