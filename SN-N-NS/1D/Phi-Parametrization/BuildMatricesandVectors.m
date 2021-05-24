%Procedure for efficient calculation of element matrices and vectors in the
%Usadel weak formulation
function [S1,S2,S3,f] = BuildMatricesandVectors(GI,weights,gamma,omega,phase,phi,Delta_0)
    %initialize matrices and vectors
    n = GI.ntot;
    S1 = sparse(n,n);
    S2 = sparse(n,n);
    S3 = sparse(n,n);
    f = zeros(n,1);
    
    %initialize function values at quadrature points in each element
    Phi = (phi(GI.elmat)*GI.PhiBS);
    G = omega./sqrt(omega^2+abs(Phi).^2);
    %calculate basis function data
    Gammainv = 1/gamma*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2); %1/gamma only in L region
    GradPhi = GI.PhiBgradS'+(GI.PhiBgradL'-GI.PhiBgradS').*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2); %basis function gradients
    h = GI.hS+(GI.hL-GI.hS)*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2);
    %build element vectors and matrices loop
    for ind1 = 1:GI.topology
        felem = GenerateVecElementVector(GI,weights,Gammainv,omega,G,phase,Delta_0,h,ind1);
        for ind2 = 1:GI.topology
            [S1elem,S2elem,S3elem] = GenerateVecElementMatrix(GI,weights,Gammainv,omega,G,phase,Delta_0,GradPhi,h,ind1,ind2);
            S1 = S1+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S1elem,n,n);
            S2 = S2+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S2elem,n,n);
            S3 = S3+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S3elem,n,n);
        end
        f = f+sparse(GI.elmat(:,ind1)',1,felem,n,1);
    end
end