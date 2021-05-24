function [S1elem,S2elem,S3elem] = GenerateVecElementMatrix(GI,weights,Gammainv,omega,G,phase,Delta_0,GradPhi,h,ind1,ind2)
    PhiS = [Delta_0*exp(-1i*phase/2)*ones(GI.nL-1,1);zeros(GI.nS-1,1);Delta_0*exp(1i*phase/2)*ones(GI.nL-1,1)];
    Gs = omega./sqrt(omega^2+abs(PhiS).^2);
    
    S1elem = h/2.*(G.^2.*GradPhi(:,ind1).*GradPhi(:,ind2))*weights;
    S2elem = h/2.*(omega*G.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    S3elem = h/2.*(Gs.*Gammainv.*G.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
end