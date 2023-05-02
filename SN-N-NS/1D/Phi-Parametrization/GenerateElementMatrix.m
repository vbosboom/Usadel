%Compute the matrices obtained in the weak formulation of 1D Usadel 
% equation in the phi-parametrization.
%S1 corresponds to the term d/dx(G_N^2 dPhi_N/dx)
%S2 corresponds to the term omega_n*Phi_N*G_N
%S3 corresponds to the term (G_S*G_N*Phi_N)/Gamma_BM

function [S1elem,S2elem,S3elem] = GenerateElementMatrix(GI,weights,Gammainv,omega,G,phase,Delta_0,GradPhi,h,ind1,ind2)
    PhiS = [Delta_0*exp(-1i*phase/2)*ones(GI.nL-1,1);zeros(GI.nS-1,1);Delta_0*exp(1i*phase/2)*ones(GI.nL-1,1)];
    Gs = omega./sqrt(omega^2+abs(PhiS).^2);
    
    S1elem = h/2.*(G.^2.*GradPhi(:,ind1).*GradPhi(:,ind2))*weights;
    S2elem = h/2.*(omega*G.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    S3elem = h/2.*(Gs.*Gammainv.*G.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
end