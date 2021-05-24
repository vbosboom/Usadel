function [S12elem,S22elem,JS22Telem,Jf12Telem,Jf12Celem,Jf22Telem,Jf32Telem,Jf32Celem,Jf42Telem,Jf42Celem] = GenerateVecElementMatrixS2N(GI,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S12elem = ksi^2*abs(GI.Deltas2)/2.*(GI.bs2(:,ind1).*GI.bs2(:,ind2)+GI.cs2(:,ind1).*GI.cs2(:,ind2));
    S22elem = ksi^2*abs(GI.Deltas2)/2.*(sin(Theta).^2).'.*(GI.bs2(:,ind1).*GI.bs2(:,ind2)+GI.cs2(:,ind1).*GI.cs2(:,ind2))*weights.';
    JS22Telem = ksi^2*abs(GI.Deltas2)/2.*((2*sin(Theta).*cos(Theta)).'.*((GI.bs2(:,ind1).*Chigrad1+GI.cs2(:,ind1).*Chigrad2).*Phi2))*weights.';
    Jf12Telem = -ksi^2*abs(GI.Deltas2)/2.*(cos(Theta).^2-sin(Theta).^2).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1.*Phi2*weights.';
    Jf12Celem = -ksi^2*abs(GI.Deltas2)/2.*((2*sin(Theta).*cos(Theta)).'.*(GI.bs2(:,ind2).*Chigrad1+GI.cs2(:,ind2).*Chigrad2).*Phi1)*weights.';
    Jf22Telem = abs(GI.Deltas2)/2.*(1i*E*cos(Theta).').*Phi1.*Phi2*weights.';
    Jf32Telem = -abs(GI.Deltas2)/2.*(0.5*sin(Theta).*(Delta.*exp(-1i*Chi)+conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf32Celem = abs(GI.Deltas2)/2.*(0.5*cos(Theta).*(-1i*Delta.*exp(-1i*Chi)+1i*conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf42Telem = -1i/2*abs(GI.Deltas2)/2.*(cos(Theta).*(Delta.*exp(-1i*Chi)-conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf42Celem = -1i/2*abs(GI.Deltas2)/2.*(sin(Theta).*(-1i*Delta.*exp(-1i*Chi)-1i*conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
end