function [S1melem,S1pelem,S2melem,S2pelem,S3melem,S3pelem] = GenerateElementMatrix(GI,weights,Gammainv,omega,omega_t,G_F,phase,Delta_0,GradPhi,h,ind1,ind2)

    PhiS = [Delta_0*exp(-1i*phase/2)*ones(GI.nL-1,1);zeros(GI.nS-1,1);Delta_0*exp(1i*phase/2)*ones(GI.nL-1,1)];
    G_S = omega./sqrt(omega^2+abs(PhiS).^2);
    
    S1melem = h/2.*(G_F.^2.*GradPhi(:,ind1).*GradPhi(:,ind2))*weights;
    S1pelem = h/2.*(G_F.^2.*GradPhi(:,ind1).*GradPhi(:,ind2))*weights;
    
    S2melem = h/2.*(omega_t*G_F.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    S2pelem = h/2.*(omega_t*G_F.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    
    S3melem = h/2.*(G_S.*Gammainv.*G_F.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    S3pelem = h/2.*(G_S.*Gammainv.*G_F.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
end