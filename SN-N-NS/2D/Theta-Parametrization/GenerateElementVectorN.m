%Compute the vectors obtained in the weak formulation of the 2D Usadel 
%equation in the theta,chi-parametrization in the N layer
%f1n corresponds to -*sin(theta_N)*cos(theta_N)*(nabla chi_N)^2
%f2n corresponds to iE*sin(theta_N)
function [f1nelem,f2nelem] = GenerateElementVectorN(GI,points,weights,E,Theta,Chigrad1,Chigrad2,ind1)
    Phi1 = GeneratePhi(points,ind1).';
    f1nelem = -abs(GI.Deltan)/2.*(sin(Theta).*cos(Theta)).'.*(Chigrad1.*Chigrad1+Chigrad2.*Chigrad2).*Phi1*weights.';
    f2nelem = abs(GI.Deltan)/2.*(1i*E*sin(Theta)).'.*Phi1*weights.';
end