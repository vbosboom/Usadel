function felem = GenerateVecElementVectorS1(GI,points,weights,omega,G,Delta,ind1)
    Phi1 = GeneratePhi(points,ind1)';
    felem = abs(GI.Deltas1)/2.*(omega*Delta.*G).'.*Phi1*weights.';
end