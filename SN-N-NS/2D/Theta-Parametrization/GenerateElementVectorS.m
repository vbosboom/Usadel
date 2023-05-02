%Compute the vectors obtained in the weak formulation of the 2D Usadel 
%equation in the theta,chi-parametrization in the S layer
%f1 corresponds to -ksi^2*sin(theta_S)*cos(theta_S)*(nabla chi_S)^2
%f2 corresponds to iE*sin(theta_S)
%f3 corresponds to 0.5*cos(theta_S)*(Delta*exp(-i*chi_S)+conj(Delta)*exp(i*chi_S))
%f4 corresponds to -0.5*1*sin(theta_S)*(Delta*exp(i*chi_S)-conj(Delta)*exp(i*chi_S))
function [f1elem,f2elem,f3elem,f4elem] = GenerateElementVectorS(Deltas,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1)
    Phi1 = GeneratePhi(points,ind1).';
    f1elem = -ksi^2.*abs(Deltas)/2.*(sin(Theta).*cos(Theta)).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1*weights.';
    f2elem = abs(Deltas)/2.*(1i*E*sin(Theta)).'.*Phi1*weights.';
    f3elem = abs(Deltas)/2.*(0.5*cos(Theta).*(Delta.*exp(-1i*Chi)+conj(Delta).*exp(1i*Chi))).'.*Phi1*weights.';
    f4elem = abs(Deltas)/2.*(-1i/2*sin(Theta).*(Delta.*exp(-1i*Chi)-conj(Delta).*exp(1i*Chi))).'.*Phi1*weights.';
end