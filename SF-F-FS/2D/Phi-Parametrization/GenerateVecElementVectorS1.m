function [fmelem,fpelem] = GenerateVecElementVectorS1(GI,points,weights,omega,G,Delta,ind1)
    
    Phi1 = GeneratePhi(points,ind1)';
    
    fmelem = abs(GI.Deltas1)/2.*(omega*conj(Delta).*G).'.*Phi1*weights.';
    fpelem = abs(GI.Deltas1)/2.*(omega*Delta.*G).'.*Phi1*weights.';
end