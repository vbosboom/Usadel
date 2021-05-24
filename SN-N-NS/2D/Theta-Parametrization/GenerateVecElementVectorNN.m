function [f1nelem,f2nelem] = GenerateVecElementVectorNN(GI,points,weights,E,Theta,Chi,Chigrad1,Chigrad2,ind1)
    Phi1 = GeneratePhi(points,ind1).';
    f1nelem = -abs(GI.Deltan)/2.*(sin(Theta).*cos(Theta)).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1*weights.';
    f2nelem = abs(GI.Deltan)/2.*(1i*E*sin(Theta)).'.*Phi1*weights.';
end