function [f11elem,f21elem,f31elem,f41elem] = GenerateVecElementVectorS1N(GI,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1)
    Phi1 = GeneratePhi(points,ind1).';
    f11elem = -abs(GI.Deltas1)/2.*(sin(Theta).*cos(Theta)).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1*weights.';
    f21elem = abs(GI.Deltas1)/2.*(1i*E*sin(Theta)).'.*Phi1*weights.';
    f31elem = abs(GI.Deltas1)/2.*(0.5*cos(Theta).*(Delta.*exp(-1i*Chi)+conj(Delta).*exp(1i*Chi))).'.*Phi1*weights.';
    f41elem = abs(GI.Deltas1)/2.*(-1i/2*sin(Theta).*(Delta.*exp(-1i*Chi)-conj(Delta).*exp(1i*Chi))).'.*Phi1*weights.';
end