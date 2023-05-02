%Modify the nonlinear sytem of matrices to take into account the bulk 
%boundary conditions at the edges of the S layers
function [S1tot,f1,S2tot,f2] = ProcessEssentialBoundaryConditions(GI,S1tot,f1,S2tot,f2,phase,Delta_0)
    
    %Essential boundary conditions for S1
    for i = 1:length(GI.ess_bc_index_s1)
        k = GI.convs1(GI.ess_bc_index_s1(i));
        f1 = f1- (S1tot(:,k))*Delta_0*exp(-1i*phase/2);
    end
    
    for i = 1:length(GI.ess_bc_index_s1)
        k = GI.convs1(GI.ess_bc_index_s1(i));
        S1tot(:,k) = 0;
        S1tot(k,:) = 0;
        S1tot(k,k) = 1;
        f1(k) = Delta_0*exp(-1i*phase/2);
    end
    
    %Essential boundary conditions for S2
    for i = 1:length(GI.ess_bc_index_s2)
        k = GI.convs2(GI.ess_bc_index_s2(i));
        f2 = f2- (S2tot(:,k))*Delta_0*exp(1i*phase/2);
    end

    for i = 1:length(GI.ess_bc_index_s2)
        k = GI.convs2(GI.ess_bc_index_s2(i));
        S2tot(:,k) = 0;
        S2tot(k,:) = 0;
        S2tot(k,k) = 1;
        f2(k) = Delta_0*exp(1i*phase/2);
    end
end