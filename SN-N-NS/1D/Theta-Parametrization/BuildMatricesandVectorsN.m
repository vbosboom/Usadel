%Procedure for efficient calculation of element matrices,vectors and 
%jacobians in the Usadel weak formulation
function [S1,f1,f2,f3,S2,f4,JS21,Jf11,Jf12,Jf21,Jf31,Jf32,Jf41,Jf42] = BuildMatricesandVectorsN(GI,weights,phase,E,gamma,chi,theta,Delta_0)
    
    %initialize matrices, vectors and jacobians
    n = GI.ntot;
    S1 = sparse(n,n);
    S2 = sparse(n,n);
    
    JS21 = sparse(n,n);
    Jf11 = sparse(n,n);
    Jf12 = sparse(n,n);
    Jf21 = sparse(n,n);
    Jf31 = sparse(n,n);
    Jf32 = sparse(n,n);
    Jf41 = sparse(n,n);
    Jf42 = sparse(n,n);
    
    f1 = zeros(n,1);
    f2 = zeros(n,1);
    f3 = zeros(n,1);
    f4 = zeros(n,1);
    
    
    %calculate basis function data
    Gammainv = 1/gamma*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2); %1/gamma only in L region
    GradPhi = GI.PhiBgradS'+(GI.PhiBgradL'-GI.PhiBgradS').*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2); %basis function gradients
    h = GI.hS+(GI.hL-GI.hS)*(mean(GI.x(GI.elmat),2)>GI.S/2 | mean(GI.x(GI.elmat),2)<-GI.S/2);
    %initialize function values at quadrature points in each element
    Theta = (theta(GI.elmat)*GI.PhiBS);
    Chi = (chi(GI.elmat)*GI.PhiBS);
    ChiGrad = sum(chi(GI.elmat).*GradPhi,2);    
    
    for ind1 = 1:GI.topology
        [f1elem,f2elem,f3elem,f4elem] = GenerateVecElementVector(GI,weights,E,Theta,Chi,Gammainv,ChiGrad,h,Delta_0,phase,ind1);
        for ind2 = 1:GI.topology
            [S1elem,S2elem,JS21elem,Jf11elem,Jf12elem,Jf21elem,Jf31elem,Jf32elem,Jf41elem,Jf42elem]...
                =GenerateVecElementMatrixN(GI,weights,E,Theta,Chi,Gammainv,GradPhi,ChiGrad,h,Delta_0,phase,ind1,ind2);
            S1 = S1+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S1elem,n,n);
            S2 = S2+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',S2elem,n,n);
            JS21 = JS21+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',JS21elem,n,n);
            Jf11 = Jf11+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',Jf11elem,n,n);
            Jf12 = Jf12+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',Jf12elem,n,n);
            Jf21 = Jf21+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',Jf21elem,n,n);
            Jf31 = Jf31+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',Jf31elem,n,n);
            Jf32 = Jf32+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',Jf32elem,n,n);
            Jf41 = Jf41+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',Jf41elem,n,n);
            Jf42 = Jf42+sparse(GI.elmat(:,ind1)',GI.elmat(:,ind2)',Jf42elem,n,n);
        end
        f1 = f1+sparse(GI.elmat(:,ind1)',1,f1elem,n,1);
        f2 = f2+sparse(GI.elmat(:,ind1)',1,f2elem,n,1);
        f3 = f3+sparse(GI.elmat(:,ind1)',1,f3elem,n,1);
        f4 = f4+sparse(GI.elmat(:,ind1)',1,f4elem,n,1);
    end
end