function [Jac,Stot,ftot] = ProcessEssentialBoundaryConditions(Jac,Stot,ftot,GI,phase,Delta_0,E )
    %boundary conditions for S1
    for i = 1:length(GI.ess_bc_index_s1)
        k = GI.convs1(GI.ess_bc_index_s1(i));
        ftot = ftot- (Stot(:,k))*atan(Delta_0/(-1i*E));  %for theta
        ftot = ftot- (Stot(:,k+GI.ns1))*(-phase/2); %for chi
    end

    for i = 1:length(GI.ess_bc_index_s1)
        k = GI.convs1(GI.ess_bc_index_s1(i));
        Stot(:,k) = 0;
        Stot(k,:) = 0;
        Stot(k,k) = 1;
        Jac(k,:) = 0;
        ftot(k) = atan(Delta_0/(-1i*E));

        Stot(:,k+GI.ns1) = 0;
        Stot(k+GI.ns1,:) = 0;
        Stot(k+GI.ns1,k+GI.ns1) = 1;
        Jac(k+GI.ns1,:) = 0;
        ftot(k+GI.ns1) = -phase/2;
    end

    %boundary conditions for S2
    for i = 1:length(GI.ess_bc_index_s2)
        k = 2*GI.ns1+GI.convs2(GI.ess_bc_index_s2(i));
        ftot = ftot- (Stot(:,k))*atan(Delta_0/(-1i*E));  %for theta
        ftot = ftot- (Stot(:,k+GI.ns2))*(phase/2); %for chi
    end

    for i = 1:length(GI.ess_bc_index_s2)
        k = 2*GI.ns1+GI.convs2(GI.ess_bc_index_s2(i));
        Stot(:,k) = 0;
        Stot(k,:) = 0;
        Stot(k,k) = 1;
        Jac(k,:) = 0;
        ftot(k) = atan(Delta_0/(-1i*E));

        Stot(:,k+GI.ns2) = 0;
        Stot(k+GI.ns2,:) = 0;
        Stot(k+GI.ns2,k+GI.ns2) = 1;
        Jac(k+GI.ns2,:) = 0;
        ftot(k+GI.ns2) = phase/2;
    end
end