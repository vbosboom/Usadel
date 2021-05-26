function [S1m,S1p,S2m,S2p,S3m,S3p,fm,fp] = BuildMatricesandVectors(GI,weights,omega,gamma,H,phase,phiP,phiCM,Delta_0)
    
    %initialize matrices and vectors
    n = GI.ntot;
    S1m = sparse(n,n);
    S1p = sparse(n,n);
    S2m = sparse(n,n);
    S2p = sparse(n,n);
    S3m = sparse(n,n);
    S3p = sparse(n,n);
    fm = zeros(n,1);
    fp = zeros(n,1);
    
    %initialize function values at quadrature points in each element
    PhiF_plus = (phiP(GI.elmat)*GI.PhiBS);
    PhiF_Cmin = (phiCM(GI.elmat)*GI.PhiBS);
    omega_t = omega+1i*H;
    G_F = omega_t./sqrt(omega_t^2+PhiF_plus.*PhiF_Cmin);
    
    %initialize element information
    Gammainv = 1/gamma*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2); %1/gamma only in L region
    GradPhi = GI.PhiBgradS'+(GI.PhiBgradL'-GI.PhiBgradS').*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2); %basis function gradients
    h = GI.hS+(GI.hL-GI.hS)*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2);
    %build element vectors and matrices loop
    for ind1 = 1:GI.topology
        [fmelem,fpelem] = GenerateElementVector(GI,weights,Gammainv,omega,omega_t,G_F,phase,Delta_0,h,ind1);
        for ind2 = 1:GI.topology
            [S1melem,S1pelem,S2melem,S2pelem,S3melem,S3pelem] = GenerateElementMatrix(GI,weights,Gammainv,omega,omega_t,G_F,phase,Delta_0,GradPhi,h,ind1,ind2);
            S1m = S1m+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S1melem,n,n);
            S1p = S1p+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S1pelem,n,n);
            S2m = S2m+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S2melem,n,n);
            S2p = S2p+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S2pelem,n,n);
            S3m = S3m+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S3melem,n,n);
            S3p = S3p+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S3pelem,n,n);
        end
        fm = fm+sparse(GI.elmat(:,ind1)',1,fmelem,n,1);
        fp = fp+sparse(GI.elmat(:,ind1)',1,fpelem,n,1);
    end
end