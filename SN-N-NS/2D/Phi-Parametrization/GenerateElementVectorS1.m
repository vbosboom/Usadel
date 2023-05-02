%Compute the vector obtained in the weak formulation of the 2D Usadel 
%equation in the phi-parametrization in the S1 layer
%f corresponds to omega_n*Delta*G_S
function felem = GenerateElementVectorS1(GI,points,weights,omega,G,Delta,ind1)
    Phi1 = GeneratePhi(points,ind1)';
    felem = abs(GI.Deltas1)/2.*(omega*Delta.*G).'.*Phi1*weights.';
end