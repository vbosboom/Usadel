function [S1nelem, S2nelem] = GenerateVecElementMatrixN(GI,points,weights,E,thetan,chin,Chigrad1,Chigrad2,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1).';
    Phi2 = GeneratePhi(points,ind2).';
    S1nelem = abs(GI.Deltan)/2.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2));
    S2nelem = abs(GI.Deltan)/2.*(sin(thetan).^2).'.*(GI.bn(:,ind1).*GI.bn(:,ind2)+GI.cn(:,ind1).*GI.cn(:,ind2))*weights.';
end