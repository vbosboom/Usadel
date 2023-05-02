%Procedure for calculation of element matrices, vectors and Jacobians in the
%2D Usadel weak formulation using the theta,chi-parametrization
function [S11,S21,f11,f21,f31,f41,fB11,fB21,S12,S22,f12,f22,f32,f42,fB12,fB22,S1n,S2n,f1n,f2n,fB1n,fB2n, ....
            JS21T,Jf11T,Jf11C,Jf21T,Jf31T,Jf31C,Jf41T,Jf41C,JfB11T1,JfB11C1,JfB11T2,JfB11C2,JfB21T1,JfB21C1,JfB21T2,JfB21C2, ...
            JS22T,Jf12T,Jf12C,Jf22T,Jf32T,Jf32C,Jf42T,Jf42C,JfB12T1,JfB12C1,JfB12T2,JfB12C2,JfB22T1,JfB22C1,JfB22T2,JfB22C2, ...
            JS23T,Jf13T,Jf13C,Jf23T,JfB1nT1S1,JfB1nC1S1,JfB1nT2S1,JfB1nC2S1,JfB2nT1S1,JfB2nC1S1,JfB2nT2S1,JfB2nC2S1, ...
            JfB1nT1S2,JfB1nC1S2,JfB1nT2S2,JfB1nC2S2,JfB2nT1S2,JfB2nC1S2,JfB2nT2S2,JfB2nC2S2] = BuildMatricesandVectorsN(GI,points,weights,points1D,weights1D,gamma_B,gamma,E,ksi,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2)
    
    %number of vertices
    n1 = GI.ns1;
    n2 = GI.ns2;
    nn = GI.nn;

    %initialize element matrices
    %inside S1
    S11 = sparse(n1,n1);
    S21 = sparse(n1,n1);
    JS21T = sparse(n1,n1);
    Jf11T = sparse(n1,n1);
    Jf11C = sparse(n1,n1);
    Jf21T = sparse(n1,n1);
    Jf31T = sparse(n1,n1);
    Jf31C = sparse(n1,n1);
    Jf41T = sparse(n1,n1);
    Jf41C = sparse(n1,n1);
    %inside S2
    S12 = sparse(n2,n2);
    S22 = sparse(n2,n2);
    JS22T = sparse(n2,n2);
    Jf12T = sparse(n2,n2);
    Jf12C = sparse(n2,n2);
    Jf22T = sparse(n2,n2);
    Jf32T = sparse(n2,n2);
    Jf32C = sparse(n2,n2);
    Jf42T = sparse(n2,n2);
    Jf42C = sparse(n2,n2);
    %inside N
    S1n = sparse(nn,nn);
    S2n = sparse(nn,nn);
    JS23T = sparse(nn,nn);
    Jf13T = sparse(nn,nn);
    Jf13C = sparse(nn,nn);
    Jf23T = sparse(nn,nn);
    %on S1-N boundary
    JfB11T1 = sparse(n1,n1);
    JfB11C1 = sparse(n1,n1);
    JfB11T2 = sparse(n1,nn);
    JfB11C2 = sparse(n1,nn);
    JfB21T1 = sparse(n1,n1);
    JfB21C1 = sparse(n1,n1);
    JfB21T2 = sparse(n1,nn);
    JfB21C2 = sparse(n1,nn);
    
    JfB1nT1S1 = sparse(nn,n1);
    JfB1nC1S1 = sparse(nn,n1);
    JfB1nT2S1 = sparse(nn,nn);
    JfB1nC2S1 = sparse(nn,nn);
    JfB2nT1S1 = sparse(nn,n1);
    JfB2nC1S1 = sparse(nn,n1);
    JfB2nT2S1 = sparse(nn,nn);
    JfB2nC2S1 = sparse(nn,nn);
    %on the S2-N boundary
    JfB12T1 = sparse(n2,n2);
    JfB12C1 = sparse(n2,n2);
    JfB12T2 = sparse(n2,nn);
    JfB12C2 = sparse(n2,nn);
    JfB22T1 = sparse(n2,n2);
    JfB22C1 = sparse(n2,n2);
    JfB22T2 = sparse(n2,nn);
    JfB22C2 = sparse(n2,nn);
    
    JfB1nT1S2 = sparse(nn,n2);
    JfB1nC1S2 = sparse(nn,n2);
    JfB1nT2S2 = sparse(nn,nn);
    JfB1nC2S2 = sparse(nn,nn);
    JfB2nT1S2 = sparse(nn,n2);
    JfB2nC1S2 = sparse(nn,n2);
    JfB2nT2S2 = sparse(nn,nn);
    JfB2nC2S2 = sparse(nn,nn);

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
    
    
    %function interpolations in S1 layer
    Theta = points*(thetas1(GI.convs1(GI.elmats1))).';
    Chi = points*(chis1(GI.convs1(GI.elmats1))).';
    Delta = points*(deltas1(GI.convs1(GI.elmats1))).';

    %Gradient interpolations in S1 layer
    Chigrad1 = sum(chis1(GI.convs1(GI.elmats1)).*GI.bs1,2);
    Chigrad2 = sum(chis1(GI.convs1(GI.elmats1)).*GI.cs1,2);

    %Calculate elements and matrices in the S1 layer
    for ind1 = 1:GI.topology
        [f11elem,f21elem,f31elem,f41elem] = GenerateElementVectorS(GI.Deltas1,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1);
        for ind2 = 1:GI.topology
            [S11elem,S21elem,JS21Telem,Jf11Telem,Jf11Celem,Jf21Telem,Jf31Telem,Jf31Celem,Jf41Telem,Jf41Celem] = GenerateElementMatrixSN(GI.Deltas1,GI.bs1,GI.cs1,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1,ind2);
            S11 = S11+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',S11elem,n1,n1);
            S21 = S21+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',S21elem,n1,n1);
            JS21T = JS21T+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',JS21Telem,n1,n1);
            Jf11T = Jf11T+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',Jf11Telem,n1,n1);
            Jf11C = Jf11C+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',Jf11Celem,n1,n1);
            Jf21T = Jf21T+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',Jf21Telem,n1,n1);
            Jf31T = Jf31T+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',Jf31Telem,n1,n1);
            Jf31C = Jf31C+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',Jf31Celem,n1,n1);
            Jf41T = Jf41T+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',Jf41Telem,n1,n1);
            Jf41C = Jf41C+sparse(GI.convs1(GI.elmats1(:,ind1))',GI.convs1(GI.elmats1(:,ind2))',Jf41Celem,n1,n1);
        end
        f11 = f11+sparse(GI.convs1(GI.elmats1(:,ind1))',1,f11elem,n1,1);
        f21 = f21+sparse(GI.convs1(GI.elmats1(:,ind1))',1,f21elem,n1,1);
        f31 = f31+sparse(GI.convs1(GI.elmats1(:,ind1))',1,f31elem,n1,1);
        f41 = f41+sparse(GI.convs1(GI.elmats1(:,ind1))',1,f41elem,n1,1);
    end
    
    
    %function interpolations in S2 layer
    Theta = points*(thetas2(GI.convs2(GI.elmats2))).';
    Chi = points*(chis2(GI.convs2(GI.elmats2))).';
    Delta = points*(deltas2(GI.convs2(GI.elmats2))).';
    
    %Gradient interpolations in S2 layer
    Chigrad1 = sum(chis2(GI.convs2(GI.elmats2)).*GI.bs2,2);
    Chigrad2 = sum(chis2(GI.convs2(GI.elmats2)).*GI.cs2,2);
    
    %Calculate elements and matrices in the S2 layer
    for ind1 = 1:GI.topology
        [f12elem,f22elem,f32elem,f42elem] = GenerateElementVectorS(GI.Deltas2,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1);
        for ind2 = 1:GI.topology
            [S12elem,S22elem,JS22Telem,Jf12Telem,Jf12Celem,Jf22Telem,Jf32Telem,Jf32Celem,Jf42Telem,Jf42Celem] = GenerateElementMatrixSN(GI.Deltas2,GI.bs2,GI.cs2,points,weights,E,ksi,Theta,Chi,Delta,Chigrad1,Chigrad2,ind1,ind2);
            S12 = S12+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',S12elem,n2,n2);
            S22 = S22+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',S22elem,n2,n2);
            JS22T = JS22T+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',JS22Telem,n2,n2);
            Jf12T = Jf12T+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',Jf12Telem,n2,n2);
            Jf12C = Jf12C+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',Jf12Celem,n2,n2);
            Jf22T = Jf22T+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',Jf22Telem,n2,n2);
            Jf32T = Jf32T+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',Jf32Telem,n2,n2);
            Jf32C = Jf32C+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',Jf32Celem,n2,n2);
            Jf42T = Jf42T+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',Jf42Telem,n2,n2);
            Jf42C = Jf42C+sparse(GI.convs2(GI.elmats2(:,ind1))',GI.convs2(GI.elmats2(:,ind2))',Jf42Celem,n2,n2);
        end
        f12 = f12+sparse(GI.convs2(GI.elmats2(:,ind1))',1,f12elem,n2,1);
        f22 = f22+sparse(GI.convs2(GI.elmats2(:,ind1))',1,f22elem,n2,1);
        f32 = f32+sparse(GI.convs2(GI.elmats2(:,ind1))',1,f32elem,n2,1);
        f42 = f42+sparse(GI.convs2(GI.elmats2(:,ind1))',1,f42elem,n2,1);
    end
    

   
    %function interpolations in N layer
    Theta = points*(thetan(GI.convn(GI.elmatn))).';

    %Gradient interpolation in N layer
    Chigrad1 = sum(chin(GI.convn(GI.elmatn)).*GI.bn,2);
    Chigrad2 = sum(chin(GI.convn(GI.elmatn)).*GI.cn,2);
    
    %Calculate elements and matrices in the N layer
    for ind1 = 1:GI.topology
        [f1nelem,f2nelem] = GenerateElementVectorN(GI,points,weights,E,Theta,Chigrad1,Chigrad2,ind1);
        for ind2 = 1:GI.topology
            [S1nelem, S2nelem, JS23Telem,Jf13Telem,Jf13Celem,Jf23Telem] = GenerateElementMatrixNN(GI,points,weights,E,Theta,Chigrad1,Chigrad2,ind1,ind2);
                S1n = S1n+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',S1nelem,nn,nn);
                S2n = S2n+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',S2nelem,nn,nn);
                JS23T = JS23T+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',JS23Telem,nn,nn);
                Jf13T = Jf13T+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',Jf13Telem,nn,nn);
                Jf13C = Jf13C+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',Jf13Celem,nn,nn);
                Jf23T = Jf23T+sparse(GI.convn(GI.elmatn(:,ind1))',GI.convn(GI.elmatn(:,ind2))',Jf23Telem,nn,nn);
        end
        f1n = f1n+sparse(GI.convn(GI.elmatn(:,ind1))',1,f1nelem,nn,1);
        f2n = f2n+sparse(GI.convn(GI.elmatn(:,ind1))',1,f2nelem,nn,1);
    end
    
    %function interpolation on S1-N boundary
    ThetaN = points1D*(thetan(GI.convn(GI.elmatbnds1))).';
    ThetaS = points1D*(thetas1(GI.convs1(GI.elmatbnds1))).';
    ChiN = points1D*(chin(GI.convn(GI.elmatbnds1))).';
    ChiS = points1D*(chis1(GI.convs1(GI.elmatbnds1))).';
    lek = sqrt((GI.x(GI.elmatbnds1(:,2))-GI.x(GI.elmatbnds1(:,1))).^2+(GI.y(GI.elmatbnds1(:,2))-GI.y(GI.elmatbnds1(:,1))).^2).';

    %calculate element matrices and vectors on S1-N boundary
    for ind1 = 1:GI.topologybnd
        [fB11elem,fB21elem,fB1nelem,fB2nelem] = GenerateBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,ThetaS,ChiS,ThetaN,ChiN,ind1);
        for ind2 = 1:GI.topologybnd
            [JfB11T1elem,JfB11C1elem,JfB11T2elem,JfB11C2elem,JfB21T1elem,JfB21C1elem,JfB21T2elem,JfB21C2elem, ...
            JfB1nT1elem,JfB1nC1elem,JfB1nT2elem,JfB1nC2elem,JfB2nT1elem,JfB2nC1elem,JfB2nT2elem,JfB2nC2elem] = GenerateBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,ThetaS,ChiS,ThetaN,ChiN,ind1,ind2);
            
            JfB11T1 = JfB11T1 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',JfB11T1elem,n1,n1);
            JfB11C1 = JfB11C1 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',JfB11C1elem,n1,n1);
            JfB11T2 = JfB11T2 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',JfB11T2elem,n1,nn);
            JfB11C2 = JfB11C2 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',JfB11C2elem,n1,nn);
            JfB21T1 = JfB21T1 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',JfB21T1elem,n1,n1);
            JfB21C1 = JfB21C1 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',JfB21C1elem,n1,n1);
            JfB21T2 = JfB21T2 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',JfB21T2elem,n1,nn);
            JfB21C2 = JfB21C2 + sparse(GI.convs1(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',JfB21C2elem,n1,nn);
            
            JfB1nT1S1 = JfB1nT1S1 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',JfB1nT1elem,nn,n1);
            JfB1nC1S1 = JfB1nC1S1 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',JfB1nC1elem,nn,n1);
            JfB1nT2S1 = JfB1nT2S1 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',JfB1nT2elem,nn,nn);
            JfB1nC2S1 = JfB1nC2S1 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',JfB1nC2elem,nn,nn);
            JfB2nT1S1 = JfB2nT1S1 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',JfB2nT1elem,nn,n1);
            JfB2nC1S1 = JfB2nC1S1 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convs1(GI.elmatbnds1(:,ind2))',JfB2nC1elem,nn,n1);
            JfB2nT2S1 = JfB2nT2S1 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',JfB2nT2elem,nn,nn);
            JfB2nC2S1 = JfB2nC2S1 + sparse(GI.convn(GI.elmatbnds1(:,ind1))',GI.convn(GI.elmatbnds1(:,ind2))',JfB2nC2elem,nn,nn);
        end
        fB11 = fB11+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',1,fB11elem,n1,1);
        fB21 = fB21+sparse(GI.convs1(GI.elmatbnds1(:,ind1))',1,fB21elem,n1,1);
        fB1n = fB1n+sparse(GI.convn(GI.elmatbnds1(:,ind1))',1,fB1nelem,nn,1);
        fB2n = fB2n+sparse(GI.convn(GI.elmatbnds1(:,ind1))',1,fB2nelem,nn,1);
    end
    
    
    
    %function interpolation on S2-N boundary
    ThetaN = points1D*(thetan(GI.convn(GI.elmatbnds2))).';
    ThetaS = points1D*(thetas2(GI.convs2(GI.elmatbnds2))).';
    ChiN = points1D*(chin(GI.convn(GI.elmatbnds2))).';
    ChiS = points1D*(chis2(GI.convs2(GI.elmatbnds2))).';
    lek = sqrt((GI.x(GI.elmatbnds2(:,2))-GI.x(GI.elmatbnds2(:,1))).^2+(GI.y(GI.elmatbnds2(:,2))-GI.y(GI.elmatbnds2(:,1))).^2).';

    %calculate element matrices and vectors on S2-N boundary
    for ind1 = 1:GI.topologybnd
        [fB12elem,fB22elem,fB1nelem,fB2nelem] = GenerateBoundaryElementVectorSN(points1D,weights1D,gamma_B,gamma,ksi,lek,ThetaS,ChiS,ThetaN,ChiN,ind1);
        for ind2 = 1:GI.topologybnd
            [JfB12T1elem,JfB12C1elem,JfB12T2elem,JfB12C2elem,JfB22T1elem,JfB22C1elem,JfB22T2elem,JfB22C2elem, ...
            JfB1nT1elem,JfB1nC1elem,JfB1nT2elem,JfB1nC2elem,JfB2nT1elem,JfB2nC1elem,JfB2nT2elem,JfB2nC2elem] = GenerateBoundaryElementMatrixSN(points1D,weights1D,gamma_B,gamma,ksi,lek,ThetaS,ChiS,ThetaN,ChiN,ind1,ind2);
            
            JfB12T1 = JfB12T1 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',JfB12T1elem,n2,n2);
            JfB12C1 = JfB12C1 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',JfB12C1elem,n2,n2);
            JfB12T2 = JfB12T2 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',JfB12T2elem,n2,nn);
            JfB12C2 = JfB12C2 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',JfB12C2elem,n2,nn);
            JfB22T1 = JfB22T1 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',JfB22T1elem,n2,n2);
            JfB22C1 = JfB22C1 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',JfB22C1elem,n2,n2);
            JfB22T2 = JfB22T2 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',JfB22T2elem,n2,nn);
            JfB22C2 = JfB22C2 + sparse(GI.convs2(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',JfB22C2elem,n2,nn);
            
            JfB1nT1S2 = JfB1nT1S2 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',JfB1nT1elem,nn,n2);
            JfB1nC1S2 = JfB1nC1S2 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',JfB1nC1elem,nn,n2);
            JfB1nT2S2 = JfB1nT2S2 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',JfB1nT2elem,nn,nn);
            JfB1nC2S2 = JfB1nC2S2 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',JfB1nC2elem,nn,nn);
            JfB2nT1S2 = JfB2nT1S2 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',JfB2nT1elem,nn,n2);
            JfB2nC1S2 = JfB2nC1S2 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convs2(GI.elmatbnds2(:,ind2))',JfB2nC1elem,nn,n2);
            JfB2nT2S2 = JfB2nT2S2 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',JfB2nT2elem,nn,nn);
            JfB2nC2S2 = JfB2nC2S2 + sparse(GI.convn(GI.elmatbnds2(:,ind1))',GI.convn(GI.elmatbnds2(:,ind2))',JfB2nC2elem,nn,nn);
        end
        fB12 = fB12+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',1,fB12elem,n2,1);
        fB22 = fB22+sparse(GI.convs2(GI.elmatbnds2(:,ind1))',1,fB22elem,n2,1);
        fB1n = fB1n+sparse(GI.convn(GI.elmatbnds2(:,ind1))',1,fB1nelem,nn,1);
        fB2n = fB2n+sparse(GI.convn(GI.elmatbnds2(:,ind1))',1,fB2nelem,nn,1);
    end
end