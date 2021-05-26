function [f1mBelem,f1pBelem,f3mBelem,f3pBelem] = GenerateVecBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,PhiS_CM,PhiS_P,G_S,PhiF_CM,PhiF_P,G_F,omega,omega_T,ind1)

    Phi1 = GeneratePhi(points1D,ind1)';
    
    f1mBelem = ksi*gamma/gamma_B*lek/2*omega/omega_T.*(G_S.*G_F.*PhiF_CM).'.*Phi1*weights1D.';
    f1pBelem = ksi*gamma/gamma_B*lek/2*omega/omega_T.*(G_S.*G_F.*PhiF_P).'.*Phi1*weights1D.';
    
    f3mBelem = 1/gamma_B*lek/2*omega_T/omega.*(G_F.*G_S.*PhiS_CM).'.*Phi1*weights1D.';
    f3pBelem = 1/gamma_B*lek/2*omega_T/omega.*(G_F.*G_S.*PhiS_P).'.*Phi1*weights1D.';
end