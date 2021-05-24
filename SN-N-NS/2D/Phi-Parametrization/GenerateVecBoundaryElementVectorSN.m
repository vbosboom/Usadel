function [f1Belem,f3Belem] = GenerateVecBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,PhiS,Gs,PhiN,Gn,ind1)
    Phi1 = GeneratePhi(points1D,ind1)';
    f1Belem = ksi*gamma/gamma_B*lek/2.*(Gs.*Gn.*PhiN).'.*Phi1*weights1D.';
    f3Belem = 1/gamma_B*lek/2.*(Gn.*Gs.*PhiS).'.*Phi1*weights1D.';
end