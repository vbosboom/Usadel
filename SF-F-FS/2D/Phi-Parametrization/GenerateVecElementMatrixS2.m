function [S1melem,S1pelem,M1melem,M1pelem] = GenerateVecElementMatrixS2(GI,points,weights,omega,ksi,G,ind1,ind2)

    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S1melem = ksi^2*abs(GI.Deltas2)/2.*(G.^2).'.*(GI.bs2(:,ind1).*GI.bs2(:,ind2)+GI.cs2(:,ind1).*GI.cs2(:,ind2))*weights.';
    S1pelem = ksi^2*abs(GI.Deltas2)/2.*(G.^2).'.*(GI.bs2(:,ind1).*GI.bs2(:,ind2)+GI.cs2(:,ind1).*GI.cs2(:,ind2))*weights.';
    M1melem = abs(GI.Deltas2)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
    M1pelem = abs(GI.Deltas2)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
end