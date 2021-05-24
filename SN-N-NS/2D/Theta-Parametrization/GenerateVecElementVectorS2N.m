function [f12elem,f22elem,f32elem,f42elem] = GenerateVecElementVectorS2N(GI,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1)
    Phi1 = GeneratePhi(points,ind1).';
    f12elem = -abs(GI.Deltas2)/2.*(sin(Theta).*cos(Theta)).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1*weights.';
    f22elem = abs(GI.Deltas2)/2.*(1i*E*sin(Theta)).'.*Phi1*weights.';
    f32elem = abs(GI.Deltas2)/2.*(0.5*cos(Theta).*(Delta.*exp(-1i*Chi)+conj(Delta).*exp(1i*Chi))).'.*Phi1*weights.';
    f42elem = abs(GI.Deltas2)/2.*(-1i/2*sin(Theta).*(Delta.*exp(-1i*Chi)-conj(Delta).*exp(1i*Chi))).'.*Phi1*weights.';
end