%Compute the boundary vectors obtained in the weak formulation of the 2D Usadel 
%equation in the theta,chi-parametrization on the SN boundary
%fB11 corresponds to ksi*gamma/gamma_B*(sin(theta_N)*cos(theta_S)*cos(chi_N-chi_S)-cos(theta_N)*sin(theta_S))
%fB21 corresponds to ksi*gamma/gamma_B*(sin(theta_N)*sin(theta_S)*sin(chi_N-chi_S))
%fB1n corresponds to 1/gamma_B*(sin(theta_S)*cos(theta_N)*cos(chi_S-chi_N)-cos(theta_S)*sin(theta_N))
%fB2n corresponds to 1/gamma_B*(sin(theta_S)*sin(theta_N)*sin(chi_S-chi_N))
function [fB11elem,fB21elem,fB1nelem,fB2nelem] = GenerateBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,ThetaS,ChiS,ThetaN,ChiN,ind1)
    Phi1 = GeneratePhi(points1D,ind1)';
    
    fB11elem = ksi*gamma/gamma_B*lek/2.*(sin(ThetaN).*cos(ThetaS).*cos(ChiN-ChiS)-cos(ThetaN).*sin(ThetaS)).'.*Phi1*weights1D.';
    fB21elem = ksi*gamma/gamma_B*lek/2.*(sin(ThetaN).*sin(ThetaS).*sin(ChiN-ChiS)).'.*Phi1*weights1D.';
    fB1nelem = 1/gamma_B*lek/2.*(sin(ThetaS).*cos(ThetaN).*cos(ChiS-ChiN)-cos(ThetaS).*sin(ThetaN)).'.*Phi1*weights1D.';
    fB2nelem = 1/gamma_B*lek/2.*(sin(ThetaS).*sin(ThetaN).*sin(ChiS-ChiN)).'.*Phi1*weights1D.';
end