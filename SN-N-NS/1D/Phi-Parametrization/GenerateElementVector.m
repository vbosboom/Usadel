%Compute the vector obtained in the weak formulation of the 1D Usadel 
%equation in the phi-parametrization corresponding to the term
%(G_N*G_S*Phi_S)/Gamma_BM
function felem = GenerateElementVector(GI,weights,Gammainv,omega,G,phase,Delta_0,h,ind1)

    %Construct Greens functions in the superconducting layer
    PhiS = [Delta_0*exp(-1i*phase/2)*ones(GI.nL-1,1);zeros(GI.nS-1,1);Delta_0*exp(1i*phase/2)*ones(GI.nL-1,1)];
    Gs = omega./sqrt(omega^2+abs(PhiS).^2);

    %Compute the vector by quadrature
    felem = h/2.*(Gs.*PhiS.*Gammainv.*G.*GI.PhiBS(ind1,:))*weights;
end