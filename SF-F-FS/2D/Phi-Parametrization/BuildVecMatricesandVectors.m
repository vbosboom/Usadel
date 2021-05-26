function [S1m,S1p,M1m,M1p,f1m,f1p,S2m,S2p,M2m,M2p,f2m,f2p,S3m,S3p,M3m,M3p,f3m,f3p] = BuildVecMatricesandVectors(GI,points,weights,points1D,weights1D,gamma_B,gamma,omega,ksi,phis1P,phis1CM,deltas1,phis2P,phis2CM,deltas2,phifP,phifCM,H)

    %numbers of nodes
    n1 = GI.ns1;
    n2 = GI.ns2;
    nn = GI.nn;
    
    %declare element matices
    %in S1
    S1m = sparse(n1,n1);
    S1p = sparse(n1,n1);
    M1m = sparse(n1,n1);
    M1p = sparse(n1,n1);
    %in S2
    S2m = sparse(n2,n2);
    S2p = sparse(n2,n2);
    M2m = sparse(n2,n2);
    M2p = sparse(n2,n2);
    %in N
    S3m = sparse(nn,nn);
    S3p = sparse(nn,nn);
    M3m = sparse(nn,nn);
    M3p = sparse(nn,nn);
    %declare element vectors
    %in S1
    f1m = zeros(n1,1);
    f1p = zeros(n1,1);
    %in S2
    f2m = zeros(n2,1);
    f2p = zeros(n2,1);
    %in N
    f3m = zeros(nn,1);
    f3p = zeros(nn,1);
    
    %Interpolation of functions in S1
    Phi_P = points*(phis1P(GI.convs1(GI.elmats1))).';
    Phi_CM = points*(phis1CM(GI.convs1(GI.elmats1))).';
    G = omega./sqrt(omega^2+Phi_P.*Phi_CM);
    Delta = points*(deltas1(GI.convs1(GI.elmats1))).';
    
    for ind1 = 1:GI.topology
        [fmelem,fpelem] = GenerateVecElementVectorS1(GI,points,weights,omega,G,Delta,ind1);
        for ind2 = 1:GI.topology
            [S1melem,S1pelem,M1melem,M1pelem] = GenerateVecElementMatrixS1(GI,points,weights,omega,ksi,G,ind1,ind2);
            S1m = S1m+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',S1melem,n1,n1);
            S1p = S1p+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',S1pelem,n1,n1);
            M1m = M1m+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',M1melem,n1,n1);
            M1p = M1p+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',M1pelem,n1,n1);
        end
        f1m = f1m+sparse(GI.convs1(GI.elmats1(:,ind1))',1,fmelem,n1,1);
        f1p = f1p+sparse(GI.convs1(GI.elmats1(:,ind1))',1,fpelem,n1,1);
    end
    
    %Interpolation of functions in S2
    Phi_P = points*(phis2P(GI.convs2(GI.elmats2))).';
    Phi_CM = points*(phis2CM(GI.convs2(GI.elmats2))).';
    G = omega./sqrt(omega^2+Phi_P.*Phi_CM);
    Delta = points*(deltas2(GI.convs2(GI.elmats2))).';
    
    for ind1 = 1:GI.topology
        [fmelem,fpelem] = GenerateVecElementVectorS2(GI,points,weights,omega,G,Delta,ind1);
        for ind2 = 1:GI.topology
            [S2melem,S2pelem,M2melem,M2pelem] = GenerateVecElementMatrixS2(GI,points,weights,omega,ksi,G,ind1,ind2);
            S2m = S2m+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',S2melem,n2,n2);
            S2p = S2p+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',S2pelem,n2,n2);
            M2m = M2m+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',M2melem,n2,n2);
            M2p = M2p+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',M2pelem,n2,n2);
        end
        f2m = f2m+sparse(GI.convs2(GI.elmats2(:,ind1))',1,fmelem,n2,1);
        f2p = f2p+sparse(GI.convs2(GI.elmats2(:,ind1))',1,fpelem,n2,1);
    end
    
    %Interpolation of funtions in F
    Phi_P = points*(phifP(GI.convn(GI.elmatn))).';
    Phi_CM = points*(phifCM(GI.convn(GI.elmatn))).';
    omega_T = omega+1i*H;
    G = omega_T./sqrt(omega_T^2+Phi_P.*Phi_CM);
    
    %Calculate element matrices in F
    for ind1 = 1:GI.topology
        for ind2 = 1:GI.topology
            [S3melem,S3pelem,M3melem,M3pelem] = GenerateVecElementMatrixF(GI,points,weights,omega_T,G,ind1,ind2);
            S3m = S3m+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',S3melem,nn,nn);
            S3p = S3p+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',S3pelem,nn,nn);
            M3m = M3m+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',M3melem,nn,nn);
            M3p = M3p+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',M3pelem,nn,nn);
        end
    end
    
    %Interpolation of functions on S1-F boundary
    PhiS_CM = points1D*(phis1CM(GI.convs1(GI.elmatbnds1))).';
    PhiS_P = points1D*(phis1P(GI.convs1(GI.elmatbnds1))).';
    G_S = omega./sqrt(omega^2+PhiS_P.*PhiS_CM);
    PhiF_CM = points1D*(phifCM(GI.convn(GI.elmatbnds1))).';
    PhiF_P = points1D*(phifP(GI.convn(GI.elmatbnds1))).';
    omega_T = omega+1i*H;
    G_F = omega_T./sqrt(omega_T^2+PhiF_P.*PhiF_CM);
    
    %element lengths at S1-F boundary
    lek = sqrt((GI.x(GI.elmatbnds1(:,2))-GI.x(GI.elmatbnds1(:,1))).^2+(GI.y(GI.elmatbnds1(:,2))-GI.y(GI.elmatbnds1(:,1))).^2).';
    
    %Calculate element matrices and vectors on S1-F boundary
    for ind1=1:GI.topologybnd
        [f1mBelem,f1pBelem,f3mBelem,f3pBelem] = GenerateVecBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,PhiS_CM,PhiS_P,G_S,PhiF_CM,PhiF_P,G_F,omega,omega_T,ind1);
        for ind2 = 1:GI.topologybnd
            [M1mBelem,M1pBelem,M3mBelem,M3pBelem] = GenerateVecBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,G_S,G_F,ind1,ind2);
            M1m = M1m+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',M1mBelem,n1,n1);
            M1p = M1p+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',M1pBelem,n1,n1);
            M3m = M3m+sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',M3mBelem,nn,nn);
            M3p = M3p+sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',M3pBelem,nn,nn);
        end
        f1m = f1m+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',1,f1mBelem,n1,1);
        f1p = f1p+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',1,f1pBelem,n1,1);
        f3m = f3m+sparse(GI.convn(GI.elmatbnds1(:,ind1))',1,f3mBelem,nn,1);
        f3p = f3p+sparse(GI.convn(GI.elmatbnds1(:,ind1))',1,f3pBelem,nn,1);
    end
    
    %Interpolation of functions on S2-F boundary
    PhiS_CM = points1D*(phis2CM(GI.convs2(GI.elmatbnds2))).';
    PhiS_P = points1D*(phis2P(GI.convs2(GI.elmatbnds2))).';
    G_S = omega./sqrt(omega^2+PhiS_P.*PhiS_CM);
    PhiF_CM = points1D*(phifCM(GI.convn(GI.elmatbnds2))).';
    PhiF_P = points1D*(phifP(GI.convn(GI.elmatbnds2))).';
    omega_T = omega+1i*H;
    G_F = omega_T./sqrt(omega_T^2+PhiF_P.*PhiF_CM);
    
    %element lengths at S2-F boundary
    lek = sqrt((GI.x(GI.elmatbnds2(:,2))-GI.x(GI.elmatbnds2(:,1))).^2+(GI.y(GI.elmatbnds2(:,2))-GI.y(GI.elmatbnds2(:,1))).^2).';
    
    %Calculate element matrices and vectors on S1-F boundary
    for ind1=1:GI.topologybnd
        [f2mBelem,f2pBelem,f3mBelem,f3pBelem] = GenerateVecBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,PhiS_CM,PhiS_P,G_S,PhiF_CM,PhiF_P,G_F,omega,omega_T,ind1);
        for ind2 = 1:GI.topologybnd
            [M2mBelem,M2pBelem,M3mBelem,M3pBelem] = GenerateVecBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,G_S,G_F,ind1,ind2);
            M2m = M2m+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',M2mBelem,n2,n2);
            M2p = M2p+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',M2pBelem,n2,n2);
            M3m = M3m+sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',M3mBelem,nn,nn);
            M3p = M3p+sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',M3pBelem,nn,nn);
        end
        f2m = f2m+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',1,f2mBelem,n2,1);
        f2p = f2p+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',1,f2pBelem,n2,1);
        f3m = f3m+sparse(GI.convn(GI.elmatbnds2(:,ind1))',1,f3mBelem,nn,1);
        f3p = f3p+sparse(GI.convn(GI.elmatbnds2(:,ind1))',1,f3pBelem,nn,1);
    end
end