function [S1melem,S1pelem,M1melem,M1pelem] = GenerateVecElementMatrixS1(GI,points,weights,omega,ksi,G,ind1,ind2)

    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S1melem = ksi^2*abs(GI.Deltas1)/2.*(G.^2).'.*(GI.bs1(:,ind1).*GI.bs1(:,ind2)+GI.cs1(:,ind1).*GI.cs1(:,ind2))*weights.';
    S1pelem = ksi^2*abs(GI.Deltas1)/2.*(G.^2).'.*(GI.bs1(:,ind1).*GI.bs1(:,ind2)+GI.cs1(:,ind1).*GI.cs1(:,ind2))*weights.';
    M1melem = abs(GI.Deltas1)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
    M1pelem = abs(GI.Deltas1)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
end