function [S1elem,S1Telem,S2elem,S2Telem,S3elem,S3Telem,S4elem,S4Telem,S5elem,S5Telem,S6elem,S6Telem] = GenerateElementMatrix(GI,weights,Gammainv,omega,H,GammaF,GammaTF,NF,GradPhi,GammaFGrad,GammaTFGrad,gammaS,gammaTS,NS,h,ind1,ind2)

    S1elem = h/2.*(GradPhi(:,ind1).*GradPhi(:,ind2))*2;
    S1Telem = h/2.*(GradPhi(:,ind1).*GradPhi(:,ind2))*2;

    S2elem = -h/2.*(2*NF.*GammaTF.*GammaFGrad.*GradPhi(:,ind2).*GI.PhiBS(ind1,:))*weights;
    S2Telem = -h/2.*(2*NF.*GammaF.*GammaTFGrad.*GradPhi(:,ind2).*GI.PhiBS(ind1,:))*weights;
    
    S3elem = h/2.*((omega+1i*H).*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    S3Telem = h/2.*((omega+1i*H).*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    
    S4elem = h/2.*(Gammainv.*NS.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    S4Telem = h/2.*(Gammainv.*NS.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    
    S5elem = h/2.*(Gammainv.*(gammaS.*gammaTS).*NS.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    S5Telem = h/2.*(Gammainv.*(gammaS.*gammaTS).*NS.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    
    %maybe factor 2 here
    S6elem = -h/2.*(Gammainv.*GammaF.*gammaTS.*NS.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
    S6Telem = -h/2.*(Gammainv.*GammaTF.*gammaS.*NS.*GI.PhiBS(ind1,:).*GI.PhiBS(ind2,:))*weights;
end