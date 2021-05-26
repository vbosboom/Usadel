function [S1,S1T,S2,S2T,S3,S3T,S4,S4T,S5,S5T,S6,S6T,f,fT,...
    JS2NN,JS2NT,JS2TN,JS2TT,JS6NN,JS6TT] = BuildMatricesandVectorsN(GI,weights,omega,gammaBM,H,phase,gammaF,gammaTF,Delta_0)

    %initialize matrices and vectors
    n = GI.ntot;
    S1 = sparse(n,n);
    S1T = sparse(n,n);
    S2 = sparse(n,n);
    S2T = sparse(n,n);
    S3 = sparse(n,n);
    S3T = sparse(n,n);
    S4 = sparse(n,n);
    S4T = sparse(n,n);
    S5 = sparse(n,n);
    S5T = sparse(n,n);
    S6 = sparse(n,n);
    S6T = sparse(n,n);
    
    f = zeros(n,1);
    fT = zeros(n,1);
    %initialize Jacobians
    JS2NN = sparse(n,n);
    JS2NT = sparse(n,n);
    JS2TN = sparse(n,n);
    JS2TT = sparse(n,n);
    
    JS6NN = sparse(n,n);
    JS6TT = sparse(n,n);

    %initialize function values at quadrature points in each element
    GammaF = (gammaF(GI.elmat)*GI.PhiBS);
    GammaTF = (gammaTF(GI.elmat)*GI.PhiBS);
    NF = 1./(1-GammaF.*GammaTF);

    Gammainv = 1/gammaBM*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2); %1/gamma only in L region
    GradPhi = GI.PhiBgradS'+(GI.PhiBgradL'-GI.PhiBgradS').*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2); %basis function gradients
    h = GI.hS+(GI.hL-GI.hS)*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2);
    
    GammaFGrad = sum(gammaF(GI.elmat).*GradPhi,2);
    GammaTFGrad = sum(gammaTF(GI.elmat).*GradPhi,2);
    
    %initialize functions in superconductor
    if(imag(omega)~=0)
        sgn = sign(imag(omega));
    else
        sgn = 1;
    end
    sgn = 1;
    Delta = [Delta_0*exp(-1i*phase/2)*ones(GI.nL-1,1);zeros(GI.nS-1,1);Delta_0*exp(1i*phase/2)*ones(GI.nL-1,1)];
    gammaS = Delta./(omega+sgn*sqrt(omega^2+abs(Delta).^2));
    gammaTS = -conj(Delta)./(omega+sgn*sqrt(omega^2+abs(Delta).^2));
    %fix problem for positive energy with the square root
    if sgn<0
        gammaS = -conj(gammaS);
        gammaTS = -conj(gammaTS);
    end
    NS = 1./(1-gammaS.*gammaTS);
    
    for ind1=1:GI.topology
        [felem,fTelem] = GenerateElementVector(GI,weights,Gammainv,gammaS,gammaTS,NS,h,ind1);
        for ind2 = 1:GI.topology
            [S1elem,S1Telem,S2elem,S2Telem,S3elem,S3Telem,S4elem,S4Telem,S5elem,S5Telem,S6elem,S6Telem,...
                JS2NNelem,JS2NTelem,JS2TNelem,JS2TTelem,JS6NNelem,JS6TTelem] = GenerateElementMatrixN(GI,weights,Gammainv,omega,H,GammaF,GammaTF,NF,GradPhi,GammaFGrad,GammaTFGrad,gammaS,gammaTS,NS,h,ind1,ind2);
            S1 = S1+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S1elem,n,n);
            S1T = S1T+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S1Telem,n,n);
            S2 = S2+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S2elem,n,n);
            S2T = S2T+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S2Telem,n,n);
            S3 = S3+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S3elem,n,n);
            S3T = S3T+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S3Telem,n,n);
            S4 = S4+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S4elem,n,n);
            S4T = S4T+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S4Telem,n,n);
            S5 = S5+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S5elem,n,n);
            S5T = S5T+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S5Telem,n,n);
            S6 = S6+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S6elem,n,n);
            S6T = S6T+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S6Telem,n,n);
            
            JS2NN = JS2NN+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',JS2NNelem,n,n);
            JS2NT = JS2NT+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',JS2NTelem,n,n);
            JS2TN = JS2TN+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',JS2TNelem,n,n);
            JS2TT = JS2TT+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',JS2TTelem,n,n);
            JS6NN = JS6NN+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',JS6NNelem,n,n);
            JS6TT = JS6TT+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',JS6TTelem,n,n);
        end
        f = f+sparse(GI.elmat(:,ind1)',1,felem,n,1);
        fT = fT+sparse(GI.elmat(:,ind1)',1,fTelem,n,1);
    end
end