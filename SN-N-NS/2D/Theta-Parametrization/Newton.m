%Solves the nonlinear matrix equation from the 2D Usadel equations in 
%the theta,chi-parametrization using the Newton method
function [thetas1,chis1,thetas2,chis2,thetan,chin,success] = Newton(GI,points,weights,points1D,weights1D,phase,E,gamma_B,gamma,ksi,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2,Delta_0,itermax,tol)
    
    %settings for the Newton method
    Diff = 10;
    iter = 0;
    success = true;

    while (Diff>tol && iter<itermax)
        
        %Calculate element matrices,vectors and Jacobians
        [S11,S21,f11,f21,f31,f41,fB11,fB21,S12,S22,f12,f22,f32,f42,fB12,fB22,S1n,S2n,f1n,f2n,fB1n,fB2n, ....
            JS21T,Jf11T,Jf11C,Jf21T,Jf31T,Jf31C,Jf41T,Jf41C,JfB11T1,JfB11C1,JfB11T2,JfB11C2,JfB21T1,JfB21C1,JfB21T2,JfB21C2, ...
            JS22T,Jf12T,Jf12C,Jf22T,Jf32T,Jf32C,Jf42T,Jf42C,JfB12T1,JfB12C1,JfB12T2,JfB12C2,JfB22T1,JfB22C1,JfB22T2,JfB22C2, ...
            JS23T,Jf13T,Jf13C,Jf23T,JfB1nT1S1,JfB1nC1S1,JfB1nT2S1,JfB1nC2S1,JfB2nT1S1,JfB2nC1S1,JfB2nT2S1,JfB2nC2S1, ...
            JfB1nT1S2,JfB1nC1S2,JfB1nT2S2,JfB1nC2S2,JfB2nT1S2,JfB2nC1S2,JfB2nT2S2,JfB2nC2S2] = BuildMatricesandVectorsN(GI,points,weights,points1D,weights1D,gamma_B,gamma,E,ksi,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2);
        
        %Aseemble the total matrices and vectors
        Stot = blkdiag(S11,S21,S12,S22,S1n,S2n);
        ftot = [f11+f21+f31+fB11;f41+fB21;f12+f22+f32+fB12;f42+fB22;f1n+f2n+fB1n;fB2n];
        
        %Assemble the total Jacobian
        Jac = [-Jf11T-Jf21T-Jf31T-JfB11T1,-Jf11C-Jf31C-JfB11C1,zeros(GI.ns1,2*GI.ns2),-JfB11T2,-JfB11C2; ...
            JS21T-Jf41T-JfB21T1,-Jf41C-JfB21C1,zeros(GI.ns1,2*GI.ns2),-JfB21T2,-JfB21C2; ...
            zeros(GI.ns2,2*GI.ns1),-Jf12T-Jf22T-Jf32T-JfB12T1,-Jf12C-Jf32C-JfB12C1,-JfB12T2,-JfB12C2; ...
            zeros(GI.ns2,2*GI.ns1),JS22T-Jf42T-JfB22T1,-Jf42C-JfB22C1,-JfB22T2,-JfB22C2; ...
            -JfB1nT1S1,-JfB1nC1S1,-JfB1nT1S2,-JfB1nC1S2,-Jf13T-Jf23T-JfB1nT2S1-JfB1nT2S2,-Jf13C-JfB1nC2S1-JfB1nC2S2; ...
            -JfB2nT1S1,-JfB2nC1S1,-JfB2nT1S2,-JfB2nC1S2,JS23T-JfB2nT2S1-JfB2nT2S2,-JfB2nC2S1-JfB2nC2S2];
        
        %apply the essential boundary conditions at S1 and S2
        [Jac,Stot,ftot] = ProcessEssentialBoundaryConditions(Jac,Stot,ftot,GI,phase,Delta_0,E);
        Jac = Jac+Stot;
        
        %cost function of the system
    	h = Stot*[thetas1;chis1;thetas2;chis2;thetan;chin]-ftot;
        
        pk = lsqminnorm(Jac,h); %descend direction of the Newton method
        %Calculate appropriate step size in the direction pk
        alpha = linesearch(0.5,10,pk,norm(h),GI,points,weights,points1D,weights1D,E,gamma_B,gamma,phase,ksi,Delta_0,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2);
        
        %calculate new solution
        thetas1 = thetas1-alpha*pk(1:GI.ns1);
        chis1 = chis1-alpha*pk(GI.ns1+1:2*GI.ns1);
        thetas2 = thetas2-alpha*pk(2*GI.ns1+1:2*GI.ns1+GI.ns2);
        chis2 = chis2-alpha*pk(2*GI.ns1+GI.ns2+1:2*(GI.ns1+GI.ns2));
        thetan = thetan-alpha*pk(2*(GI.ns1+GI.ns2)+1:2*(GI.ns1+GI.ns2)+GI.nn);
        chin = chin-alpha*pk(2*(GI.ns1+GI.ns2)+GI.nn+1:end);

        %test if norm of h is sufficiently low
        Diff = norm(h);
        iter = iter+1;
        Difference(iter) = Diff;
    end
end