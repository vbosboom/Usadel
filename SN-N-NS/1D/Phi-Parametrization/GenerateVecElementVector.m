function felem = GenerateVecElementVector(GI,weights,Gammainv,omega,G,phase,Delta_0,h,ind1)
    PhiS = [Delta_0*exp(-1i*phase/2)*ones(GI.nL-1,1);zeros(GI.nS-1,1);Delta_0*exp(1i*phase/2)*ones(GI.nL-1,1)];
    Gs = omega./sqrt(omega^2+abs(PhiS).^2);
    
    felem = h/2.*(Gs.*PhiS.*Gammainv.*G.*GI.PhiBS(ind1,:))*weights;
end