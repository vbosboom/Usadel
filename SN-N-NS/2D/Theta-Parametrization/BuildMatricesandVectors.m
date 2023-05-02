%Builds the elment matrices, vectors and Jacobians needed for the Newton
%method
function [S11,S21,f11,f21,f31,f41,fB11,fB21,S12,S22,f12,f22,f32,f42,fB12,fB22,S1n,S2n,f1n,f2n,fB1n,fB2n] = BuildMatricesandVectors(GI,points,weights,points1D,weights1D,gamma_B,gamma,E,ksi,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2)
    
    %number of vertices
    n1 = GI.ns1;
    n2 = GI.ns2;
    nn = GI.nn;
    %initialize element matrices
    %inside S1
    S11 = sparse(n1,n1);
    S21 = sparse(n1,n1);
    %inside S2
    S12 = sparse(n2,n2);
    S22 = sparse(n2,n2);
    %inside N
    S1n = sparse(nn,nn);
    S2n = sparse(nn,nn);
    %initialize element vectors
    %inside S1
    f11 = zeros(n1,1);
    f21 = zeros(n1,1);
    f31 = zeros(n1,1);
    f41 = zeros(n1,1);
    %inside S2
    f12 = zeros(n2,1);
    f22 = zeros(n2,1);
    f32 = zeros(n2,1);
    f42 = zeros(n2,1);
    %inside N
    f1n = zeros(nn,1);
    f2n = zeros(nn,1);
    % On S1-N boundary
    fB11 = zeros(n1,1);
    fB21 = zeros(n1,1);
    % On S2-N boundary
    fB12 = zeros(n2,1);
    fB22 = zeros(n2,1);
    %on whole boundary
    fB1n = zeros(nn,1);
    fB2n = zeros(nn,1);
    
    %calculate element matrices and vectors in S1 layer
    %function interpolations
    Theta = points*(thetas1(GI.convs1(GI.elmats1))).';
    Chi = points*(chis1(GI.convs1(GI.elmats1))).';
    Delta = points*(deltas1(GI.convs1(GI.elmats1))).';
    %gradients
    Chigrad1 = sum(chis1(GI.convs1(GI.elmats1)).*GI.bs1,2);
    Chigrad2 = sum(chis1(GI.convs1(GI.elmats1)).*GI.cs1,2);
    for ind1 = 1:GI.topology
        [f11elem,f21elem,f31elem,f41elem] = GenerateElementVectorS(GI.Deltas1,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1);
        for ind2 = 1:GI.topology
            [S11elem,S21elem] = GenerateElementMatrixS(GI.Deltas1,GI.bs1,GI.cs1,weights,ksi,Theta,ind1,ind2);
            S11 = S11+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',S11elem,n1,n1);
            S21 = S21+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',S21elem,n1,n1);
        end
        f11 = f11+sparse(GI.convs1(GI.elmats1(:,ind1))',1,f11elem,n1,1);
        f21 = f21+sparse(GI.convs1(GI.elmats1(:,ind1))',1,f21elem,n1,1);
        f31 = f31+sparse(GI.convs1(GI.elmats1(:,ind1))',1,f31elem,n1,1);
        f41 = f41+sparse(GI.convs1(GI.elmats1(:,ind1))',1,f41elem,n1,1);
    end
    
    %calculate element matrices and vectors in S2 layer
    %function interpolations
    Theta = points*(thetas2(GI.convs2(GI.elmats2))).';
    Chi = points*(chis2(GI.convs2(GI.elmats2))).';
    Delta = points*(deltas2(GI.convs2(GI.elmats2))).';
    %gradients
    Chigrad1 = sum(chis2(GI.convs2(GI.elmats2)).*GI.bs2,2);
    Chigrad2 = sum(chis2(GI.convs2(GI.elmats2)).*GI.cs2,2);
    
    for ind1 = 1:GI.topology
        [f12elem,f22elem,f32elem,f42elem] = GenerateElementVectorS(GI.Deltas2,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1);
        for ind2 = 1:GI.topology
            [S12elem,S22elem] = GenerateElementMatrixS(GI.Deltas2,GI.bs2,GI.cs2,weights,ksi,Theta,ind1,ind2);
            S12 = S12+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',S12elem,n2,n2);
            S22 = S22+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',S22elem,n2,n2);
        end
        f12 = f12+sparse(GI.convs2(GI.elmats2(:,ind1))',1,f12elem,n2,1);
        f22 = f22+sparse(GI.convs2(GI.elmats2(:,ind1))',1,f22elem,n2,1);
        f32 = f32+sparse(GI.convs2(GI.elmats2(:,ind1))',1,f32elem,n2,1);
        f42 = f42+sparse(GI.convs2(GI.elmats2(:,ind1))',1,f42elem,n2,1);
    end
    

    %calculate element matrices and vectors in N layer
    %function interpolations
    Theta = points*(thetan(GI.convn(GI.elmatn))).';
    %gradients
    Chigrad1 = sum(chin(GI.convn(GI.elmatn)).*GI.bn,2);
    Chigrad2 = sum(chin(GI.convn(GI.elmatn)).*GI.cn,2);
    
    for ind1 = 1:GI.topology
        [f1nelem,f2nelem] = GenerateElementVectorN(GI,points,weights,E,Theta,Chigrad1,Chigrad2,ind1);
        for ind2 = 1:GI.topology
            [S1nelem, S2nelem] = GenerateElementMatrixN(GI,weights,Theta,ind1,ind2);
                S1n = S1n+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',S1nelem,nn,nn);
                S2n = S2n+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',S2nelem,nn,nn);
        end
        f1n = f1n+sparse(GI.convn(GI.elmatn(:,ind1))',1,f1nelem,nn,1);
        f2n = f2n+sparse(GI.convn(GI.elmatn(:,ind1))',1,f2nelem,nn,1);
    end
    
    %calculate element matrices and vectors on S1-N boundary
    ThetaN = points1D*(thetan(GI.convn(GI.elmatbnds1))).';
    ThetaS = points1D*(thetas1(GI.convs1(GI.elmatbnds1))).';
    ChiN = points1D*(chin(GI.convn(GI.elmatbnds1))).';
    ChiS = points1D*(chis1(GI.convs1(GI.elmatbnds1))).';
    lek = sqrt((GI.x(GI.elmatbnds1(:,2))-GI.x(GI.elmatbnds1(:,1))).^2+(GI.y(GI.elmatbnds1(:,2))-GI.y(GI.elmatbnds1(:,1))).^2).';

    for ind1 = 1:GI.topologybnd
        [fB11elem,fB21elem,fB1nelem,fB2nelem] = GenerateBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,ThetaS,ChiS,ThetaN,ChiN,ind1);
        fB11 = fB11+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',1,fB11elem,n1,1);
        fB21 = fB21+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',1,fB21elem,n1,1);
        fB1n = fB1n+sparse(GI.convn(GI.elmatbnds1(:,ind1))',1,fB1nelem,nn,1);
        fB2n = fB2n+sparse(GI.convn(GI.elmatbnds1(:,ind1))',1,fB2nelem,nn,1);
    end
    
    
    
    %calculate element matrices and vectors on S1-N boundary
    ThetaN = points1D*(thetan(GI.convn(GI.elmatbnds2))).';
    ThetaS = points1D*(thetas2(GI.convs2(GI.elmatbnds2))).';
    ChiN = points1D*(chin(GI.convn(GI.elmatbnds2))).';
    ChiS = points1D*(chis2(GI.convs2(GI.elmatbnds2))).';
    lek = sqrt((GI.x(GI.elmatbnds2(:,2))-GI.x(GI.elmatbnds2(:,1))).^2+(GI.y(GI.elmatbnds2(:,2))-GI.y(GI.elmatbnds2(:,1))).^2).';

    for ind1 = 1:GI.topologybnd
        [fB12elem,fB22elem,fB1nelem,fB2nelem] = GenerateBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,ThetaS,ChiS,ThetaN,ChiN,ind1);
        fB12 = fB12+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',1,fB12elem,n2,1);
        fB22 = fB22+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',1,fB22elem,n2,1);
        fB1n = fB1n+sparse(GI.convn(GI.elmatbnds2(:,ind1))',1,fB1nelem,nn,1);
        fB2n = fB2n+sparse(GI.convn(GI.elmatbnds2(:,ind1))',1,fB2nelem,nn,1);
    end
end