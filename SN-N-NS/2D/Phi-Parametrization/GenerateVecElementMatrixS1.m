function [S1elem,M1elem] = GenerateVecElementMatrixS1(GI,points,weights,omega,ksi,G,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S1elem = ksi^2*abs(GI.Deltas1)/2.*(G.^2).'.*(GI.bs1(:,ind1).*GI.bs1(:,ind2)+GI.cs1(:,ind1).*GI.cs1(:,ind2))*weights.';
    M1elem = abs(GI.Deltas1)/2.*(omega*G).'.*Phi1.*Phi2*weights.';
end