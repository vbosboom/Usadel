function [felem,fTelem] = GenerateElementVector(GI,weights,Gammainv,gammaS,gammaTS,NS,h,ind1)

    felem = h/2.*(Gammainv.*gammaS.*NS.*GI.PhiBS(ind1,:))*weights;
    fTelem = h/2.*(Gammainv.*gammaTS.*NS.*GI.PhiBS(ind1,:))*weights;
end