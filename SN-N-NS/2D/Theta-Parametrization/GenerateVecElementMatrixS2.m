function [S12elem,S22elem] = GenerateVecElementMatrixS2(GI,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S12elem = ksi^2*abs(GI.Deltas2)/2.*(GI.bs2(:,ind1).*GI.bs2(:,ind2)+GI.cs2(:,ind1).*GI.cs2(:,ind2));
    S22elem = ksi^2*abs(GI.Deltas2)/2.*(sin(Theta).^2).'.*(GI.bs2(:,ind1).*GI.bs2(:,ind2)+GI.cs2(:,ind1).*GI.cs2(:,ind2))*weights.';
end