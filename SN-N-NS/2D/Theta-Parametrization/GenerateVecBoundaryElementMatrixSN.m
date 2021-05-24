function [JfB11T1elem,JfB11C1elem,JfB11T2elem,JfB11C2elem,JfB21T1elem,JfB21C1elem,JfB21T2elem,JfB21C2elem, ...
            JfB1nT1elem,JfB1nC1elem,JfB1nT2elem,JfB1nC2elem,JfB2nT1elem,JfB2nC1elem,JfB2nT2elem,JfB2nC2elem] = GenerateVecBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,ThetaS,ChiS,ThetaN,ChiN,ind1,ind2)
    
    Phi1 = GeneratePhi(points1D,ind1)';
    Phi2 = GeneratePhi(points1D,ind2)';
    
    JfB11T1elem = ksi*gamma/gamma_B*lek/2.*(-sin(ThetaS).*sin(ThetaN).*cos(ChiS-ChiN)-cos(ThetaN).*cos(ThetaS)).'.*Phi1.*Phi2*weights1D.';
    JfB11C1elem = ksi*gamma/gamma_B*lek/2.*(sin(ThetaN).*cos(ThetaS).*sin(ChiN-ChiS)).'.*Phi1.*Phi2*weights1D.';
    JfB11T2elem = ksi*gamma/gamma_B*lek/2.*(cos(ThetaN).*cos(ThetaS).*cos(ChiS-ChiN)+sin(ThetaN).*sin(ThetaS)).'.*Phi1.*Phi2*weights1D.';
    JfB11C2elem = ksi*gamma/gamma_B*lek/2.*(sin(ThetaN).*cos(ThetaS).*sin(ChiS-ChiN)).'.*Phi1.*Phi2*weights1D.';
    JfB21T1elem = ksi*gamma/gamma_B*lek/2.*(sin(ThetaN).*cos(ThetaS).*sin(ChiN-ChiS)).'.*Phi1.*Phi2*weights1D.';
    JfB21C1elem = -ksi*gamma/gamma_B*lek/2.*(sin(ThetaN).*sin(ThetaS).*cos(ChiN-ChiS)).'.*Phi1.*Phi2*weights1D.';
    JfB21T2elem = ksi*gamma/gamma_B*lek/2.*(cos(ThetaN).*sin(ThetaS).*sin(ChiN-ChiS)).'.*Phi1.*Phi2*weights1D.';
    JfB21C2elem = ksi*gamma/gamma_B*lek/2.*(sin(ThetaN).*sin(ThetaS).*cos(ChiN-ChiS)).'.*Phi1.*Phi2*weights1D.';
    
    JfB1nT1elem = 1/gamma_B*lek/2.*(cos(ThetaS).*cos(ThetaN).*cos(ChiS-ChiN)+sin(ThetaS).*sin(ThetaN)).'.*Phi1.*Phi2*weights1D.';
    JfB1nC1elem = 1/gamma_B*lek/2.*(sin(ThetaS).*cos(ThetaN).*sin(ChiN-ChiS)).'.*Phi1.*Phi2*weights1D.';
    JfB1nT2elem = 1/gamma_B*lek/2.*(-sin(ThetaN).*sin(ThetaS).*cos(ChiS-ChiN)-cos(ThetaS).*cos(ThetaN)).'.*Phi1.*Phi2*weights1D.';
    JfB1nC2elem = 1/gamma_B*lek/2.*(sin(ThetaS).*cos(ThetaN).*sin(ChiS-ChiN)).'.*Phi1.*Phi2*weights1D.';
    JfB2nT1elem = 1/gamma_B*lek/2.*(sin(ThetaN).*cos(ThetaS).*sin(ChiS-ChiN)).'.*Phi1.*Phi2*weights1D.';
    JfB2nC1elem = 1/gamma_B*lek/2.*(sin(ThetaN).*sin(ThetaS).*cos(ChiS-ChiN)).'.*Phi1.*Phi2*weights1D.';
    JfB2nT2elem = 1/gamma_B*lek/2.*(cos(ThetaN).*sin(ThetaS).*sin(ChiS-ChiN)).'.*Phi1.*Phi2*weights1D.';
    JfB2nC2elem = -1/gamma_B*lek/2.*(sin(ThetaN).*sin(ThetaS).*cos(ChiS-ChiN)).'.*Phi1.*Phi2*weights1D.';
end