function [S1,M1,f1,S2,M2,f2,S3,M3,f3] = BuildVecMatricesandVectors(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1,deltas1,phis2,deltas2,phin)
    
    %numbers of nodes
    n1 = GI.ns1;
    n2 = GI.ns2;
    nn = GI.nn;
    
    %declare element matices
    %in S1
    S1 = sparse(n1,n1);
    M1 = sparse(n1,n1);
    %in S2
    S2 = sparse(n2,n2);
    M2 = sparse(n2,n2);
    %in N
    S3 = sparse(nn,nn);
    M3 = sparse(nn,nn);
    %declare element vectors
    %in S1
    f1 = zeros(n1,1);
    %in S2
    f2 = zeros(n2,1);
    %in N
    f3 = zeros(nn,1);
    
    %Interpolation of functions in S1
    Phi = points*(phis1(GI.convs1(GI.elmats1))).';
    G = omega./sqrt(omega^2+abs(Phi).^2);
    Delta = points*(deltas1(GI.convs1(GI.elmats1))).';
    
    %Calculate element matrices and vectors in S1
    for ind1 = 1:GI.topology
        felem = GenerateVecElementVectorS1(GI,points,weights,omega,G,Delta,ind1);
        for ind2 = 1:GI.topology
            [S1elem,M1elem] = GenerateVecElementMatrixS1(GI,points,weights,omega,ksi,G,ind1,ind2);
            S1 = S1+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',S1elem,n1,n1);
            M1 = M1+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',M1elem,n1,n1);
        end
        f1 = f1+sparse(GI.convs1(GI.elmats1(:,ind1))',1,felem,n1,1);
    end
    
    %Interpolation of functions in S2
    Phi = points*(phis2(GI.convs2(GI.elmats2))).';
    G = omega./sqrt(omega^2+abs(Phi).^2);
    Delta = points*(deltas2(GI.convs2(GI.elmats2))).';

    %Calculate element matrices and vectors in S2
    for ind1 = 1:GI.topology
        felem = GenerateVecElementVectorS2(GI,points,weights,omega,G,Delta,ind1);
        for ind2 = 1:GI.topology
            [S2elem,M2elem] = GenerateVecElementMatrixS2(GI,points,weights,omega,ksi,G,ind1,ind2);
            S2 = S2+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',S2elem,n2,n2);
            M2 = M2+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',M2elem,n2,n2);
        end
        f2 = f2+sparse(GI.convs2(GI.elmats2(:,ind1))',1,felem,n2,1);
    end
    
    %Interpolation of functions in N
    Phi = points*(phin(GI.convn(GI.elmatn))).';
    G = omega./sqrt(omega^2+abs(Phi).^2);
    
    %Calculate element matrices in N
    for ind1 = 1:GI.topology
        for ind2 = 1:GI.topology
            [S3elem,M3elem] = GenerateVecElementMatrixN(GI,points,weights,omega,ksi,G,ind1,ind2);
            S3 = S3+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',S3elem,nn,nn);
            M3 = M3+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',M3elem,nn,nn);
        end
    end
    
    %Interpolation of function on S1-N boundary
    PhiS = points1D*(phis1(GI.convs1(GI.elmatbnds1))).';
    Gs = omega./sqrt(omega^2+abs(PhiS).^2);
    PhiN = points1D*(phin(GI.convn(GI.elmatbnds1))).';
    Gn = omega./sqrt(omega^2+abs(PhiN).^2);
    %element lengths at S1-N boundary
    lek = sqrt((GI.x(GI.elmatbnds1(:,2))-GI.x(GI.elmatbnds1(:,1))).^2+(GI.y(GI.elmatbnds1(:,2))-GI.y(GI.elmatbnds1(:,1))).^2).';

    %Generate element matrices and vectors on S1-N boundary
    for ind1 = 1:GI.topologybnd
        [f1Belem,f3Belem] = GenerateVecBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,PhiS,Gs,PhiN,Gn,ind1);
        for ind2 = 1:GI.topologybnd
            [M1Belem,M3Belem] = GenerateVecBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,Gs,Gn,ind1,ind2);
            M1 = M1 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',M1Belem,n1,n1);
            M3 = M3 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',M3Belem,nn,nn);
        end
        f1 = f1+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',1,f1Belem,n1,1);
        f3 = f3+sparse(GI.convn(GI.elmatbnds1(:,ind1))',1,f3Belem,nn,1);
    end
    
    %Interpolation of function on S2-N boundary
    PhiS = points1D*(phis2(GI.convs2(GI.elmatbnds2))).';
    Gs = omega./sqrt(omega^2+abs(PhiS).^2);
    PhiN = points1D*(phin(GI.convn(GI.elmatbnds2))).';
    Gn = omega./sqrt(omega^2+abs(PhiN).^2);
    %Lenght of elements on S2-N boundary
    lek = sqrt((GI.x(GI.elmatbnds2(:,2))-GI.x(GI.elmatbnds2(:,1))).^2+(GI.y(GI.elmatbnds2(:,2))-GI.y(GI.elmatbnds2(:,1))).^2).';
    
    for ind1 = 1:GI.topologybnd
        [f2Belem,f3Belem] = GenerateVecBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,PhiS,Gs,PhiN,Gn,ind1);
        for ind2 = 1:GI.topologybnd
            [M2Belem,M3Belem] = GenerateVecBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,Gs,Gn,ind1,ind2);
            M2 = M2 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',M2Belem,n2,n2);
            M3 = M3 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',M3Belem,nn,nn);
        end
        f2 = f2+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',1,f2Belem,n2,1);
        f3 = f3+sparse(GI.convn(GI.elmatbnds2(:,ind1))',1,f3Belem,nn,1);
    end
end