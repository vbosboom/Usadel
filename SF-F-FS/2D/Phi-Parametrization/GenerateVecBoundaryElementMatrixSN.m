function [M1mBelem,M1pBelem,M3mBelem,M3pBelem] = GenerateVecBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,G_S,G_F,ind1,ind2)

    Phi1 = GeneratePhi(points1D,ind1)';
    Phi2 = GeneratePhi(points1D,ind2)';
    
    M1mBelem = ksi*gamma/gamma_B*lek/2.*(G_S.*G_F).'.*Phi1.*Phi2*weights1D.';
    M1pBelem = ksi*gamma/gamma_B*lek/2.*(G_S.*G_F).'.*Phi1.*Phi2*weights1D.';
    M3mBelem = 1/gamma_B*lek/2.*(G_F.*G_S).'.*Phi1.*Phi2*weights1D.';
    M3pBelem = 1/gamma_B*lek/2.*(G_F.*G_S).'.*Phi1.*Phi2*weights1D.';
end