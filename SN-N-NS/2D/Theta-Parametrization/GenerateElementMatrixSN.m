%Compute the matrices and Jacobians obtained in the weak formulation of the 2D Usadel 
%equation in the theta,chi-parametrization in the S layer
%S1 corresponds to -ksi^2*(nabla^2 theta_S)
%S2 corresponds to -ksi^2*nabla(sin(theta_S)^2*nabla chi_S)
function [S1elem,S2elem,JS2Telem,Jf1Telem,Jf1Celem,Jf2Telem,Jf3Telem,Jf3Celem,Jf4Telem,Jf4Celem] = GenerateElementMatrixSN(Deltas,b,c,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1,ind2)
    Phi1 = GeneratePhi(points,ind1)';
    Phi2 = GeneratePhi(points,ind2)';
    
    S1elem = ksi^2*abs(Deltas)/2.*(b(:,ind1).*b(:,ind2)+c(:,ind1).*c(:,ind2));
    S2elem = ksi^2*abs(Deltas)/2.*(sin(Theta).^2).'.*(b(:,ind1).*b(:,ind2)+c(:,ind1).*c(:,ind2))*weights.';
    
    JS2Telem = ksi^2*abs(Deltas)/2.*((2*sin(Theta).*cos(Theta)).'.*((b(:,ind1).*Chigrad1+c(:,ind1).*Chigrad2).*Phi2))*weights.';
    Jf1Telem = -ksi^2*abs(Deltas)/2.*(cos(Theta).^2-sin(Theta).^2).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1.*Phi2*weights.';
    Jf1Celem = -ksi^2*abs(Deltas)/2.*((2*sin(Theta).*cos(Theta)).'.*(b(:,ind2).*Chigrad1+c(:,ind2).*Chigrad2).*Phi1)*weights.';
    Jf2Telem = abs(Deltas)/2.*(1i*E*cos(Theta).').*Phi1.*Phi2*weights.';
    Jf3Telem = -abs(Deltas)/2.*(0.5*sin(Theta).*(Delta.*exp(-1i*Chi)+conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf3Celem = abs(Deltas)/2.*(0.5*cos(Theta).*(-1i*Delta.*exp(-1i*Chi)+1i*conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf4Telem = -1i/2*abs(Deltas)/2.*(cos(Theta).*(Delta.*exp(-1i*Chi)-conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
    Jf4Celem = -1i/2*abs(Deltas)/2.*(sin(Theta).*(-1i*Delta.*exp(-1i*Chi)-1i*conj(Delta).*exp(1i*Chi))).'.*Phi1.*Phi2*weights.';
end