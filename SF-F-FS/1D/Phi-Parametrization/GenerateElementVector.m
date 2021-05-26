function [fmelem,fpelem] = GenerateElementVector(GI,weights,Gammainv,omega,omega_t,G_F,phase,Delta_0,h,ind1)

    PhiS = [Delta_0*exp(-1i*phase/2)*ones(GI.nL-1,1);zeros(GI.nS-1,1);Delta_0*exp(1i*phase/2)*ones(GI.nL-1,1)];
    G_S = omega./sqrt(omega^2+abs(PhiS).^2);
    
    fmelem = h/2.*(omega_t/omega*G_S.*G_F.*conj(PhiS).*Gammainv.*GI.PhiBS(ind1,:))*weights;
    fpelem = h/2.*(omega_t/omega*G_S.*G_F.*PhiS.*Gammainv.*GI.PhiBS(ind1,:))*weights;
end