function [S11elem,S21elem,JS21Telem,Jf11Telem,Jf11Celem,Jf21Telem,Jf31Telem,Jf31Celem,Jf41Telem,Jf41Celem] = GenerateVecElementMatrixS1N(GI,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S11elem = ksi^2*abs(GI.Deltas1)/2.*(GI.bs1(:,ind1).*GI.bs1(:,ind2)+GI.cs1(:,ind1).*GI.cs1(:,ind2));
    S21elem = ksi^2*abs(GI.Deltas1)/2.*(sin(Theta).^2).'.*(GI.bs1(:,ind1).*GI.bs1(:,ind2)+GI.cs1(:,ind1).*GI.cs1(:,ind2))*weights.';
    JS21Telem = ksi^2*abs(GI.Deltas1)/2.*((2*sin(Theta).*cos(Theta)).'.*((GI.bs1(:,ind1).*Chigrad1+GI.cs1(:,ind1).*Chigrad2).*Phi2))*weights.';
    Jf11Telem = -ksi^2*abs(GI.Deltas1)/2.*(cos(Theta).^2-sin(Theta).^2).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1.*Phi2*weights.';
    Jf11Celem = -ksi^2*abs(GI.Deltas1)/2.*((2*sin(Theta).*cos(Theta)).'.*(GI.bs1(:,ind2).*Chigrad1+GI.cs1(:,ind2).*Chigrad2).*Phi1)*weights.';
    Jf21Telem = abs(GI.Deltas1)/2.*(1i*E*cos(Theta).').*Phi1.*Phi2*weights.';
    Jf31Telem = -abs(GI.Deltas1)/2.*(0.5*sin(Theta).*(Delta.*exp(-1i*Chi)+conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf31Celem = abs(GI.Deltas1)/2.*(0.5*cos(Theta).*(-1i*Delta.*exp(-1i*Chi)+1i*conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf41Telem = -1i/2*abs(GI.Deltas1)/2.*(cos(Theta).*(Delta.*exp(-1i*Chi)-conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf41Celem = -1i/2*abs(GI.Deltas1)/2.*(sin(Theta).*(-1i*Delta.*exp(-1i*Chi)-1i*conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
end