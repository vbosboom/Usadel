function [S11elem,S21elem] = GenerateVecElementMatrixS1(GI,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S11elem = ksi^2*abs(GI.Deltas1)/2.*(GI.bs1(:,ind1).*GI.bs1(:,ind2)+GI.cs1(:,ind1).*GI.cs1(:,ind2));
    S21elem = ksi^2*abs(GI.Deltas1)/2.*(sin(Theta).^2).'.*(GI.bs1(:,ind1).*GI.bs1(:,ind2)+GI.cs1(:,ind1).*GI.cs1(:,ind2))*weights.';
end