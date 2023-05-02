%Compute the boundary matrices obtained in the weak formulation of the 2D 
%Usadel equation in the phi-parametrization on the SN boundary
%M1 corresponds to ksi*gamma/gamma_B*(G_S*G_N*Phi_S)
%M3 corresponds to (G_S*G_N*Phi_N)/gamma_B
function [M1Belem,M3Belem] = GenerateBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,Gs,Gn,ind1,ind2)
    Phi1 = GeneratePhi(points1D,ind1)';
    Phi2 = GeneratePhi(points1D,ind2)';

    M1Belem = ksi*gamma/gamma_B*lek/2.*(Gs.*Gn).'.*Phi1.*Phi2*weights1D.';
    M3Belem = 1/gamma_B*lek/2.*(Gn.*Gs).'.*Phi1.*Phi2*weights1D.';
end