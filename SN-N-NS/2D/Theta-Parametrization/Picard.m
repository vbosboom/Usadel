%Solves the system of equations following from the finite element method
%using the Picard method
function [thetas1,chis1,thetas2,chis2,thetan,chin,success] = Picard(GI,points,weights,points1D,weights1D,phase,E,gamma_B,gamma,ksi,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2,Delta_0,itermax,tol)
    
    %settings for the Picard method
    Diff = 10;
    iter = 0;
    success = true;

    while (Diff>tol && iter<itermax)
        %solution from previous iteration
        thetas1old = thetas1;
        chis1old = chis1;
        thetas2old = thetas2;
        chis2old = chis2;
        thetanold = thetan;
        chinold = chin;
        
        %calculate element matrices and vectors
        [S11,S21,f11,f21,f31,f41,fB11,fB21,S12,S22,f12,f22,f32,f42,fB12,fB22,S1n,S2n,f1n,f2n,fB1n,fB2n] = BuildVecMatricesandVectors(GI,points,weights,points1D,weights1D,gamma_B,gamma,E,ksi,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2);
        Stot = blkdiag(S11,S21,S12,S22,S1n,S2n);
        ftot = [f11+f21+f31+fB11;f41+fB21;f12+f22+f32+fB12;f42+fB22;f1n+f2n+fB1n;fB2n];
        
        %apply the essential boundary conditions at S1 and S2
        [~,Stot,ftot] = ProcessEssentialBoundaryConditions(sparse(2*(GI.ns1+GI.ns2+GI.nn),2*(GI.ns1+GI.ns2+GI.nn)),Stot,ftot,GI,phase,Delta_0,E);
                
        %calculate new solution
        sol = lsqminnorm(Stot,ftot);
        thetas1 = sol(1:GI.ns1);
        chis1 = sol(GI.ns1+1:2*GI.ns1);
        thetas2 = sol(2*GI.ns1+1:2*GI.ns1+GI.ns2);
        chis2 = sol(2*GI.ns1+GI.ns2+1:2*(GI.ns1+GI.ns2));
        thetan = sol(2*(GI.ns1+GI.ns2)+1:2*(GI.ns1+GI.ns2)+GI.nn);
        chin = sol(2*(GI.ns1+GI.ns2)+GI.nn+1:end);

        %test if norm of Green's function converges
        Diff = norm(thetas1-thetas1old);
        iter = iter+1;
        Difference(iter) = Diff;
    end
end