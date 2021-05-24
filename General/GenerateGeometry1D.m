%Generate the geometry information datastructure (GI) contain information about
%the finite element mesh and basis functions

function GI = GenerateGeometry1D(S,L,nS,nL,points)
    %make the empty datastructure
    GI = struct;

    %create vector with nodal positions;
    x1 = linspace(-L/2,-S/2,nL);
    x2 = linspace(-S/2,S/2,nS);
    x3 = linspace(S/2,L/2,nL);
    x = [x1(1:end-1),x2(1:end-1),x3];
    ntot = nS+2*nL-2;
    hS = S/(nS-1);
    hL = (L/2-S/2)/(nL-1);
    
    %info about the basis functions
    GI.PhiBS = GeneratePhi1D(0,hS,hS,hS/2*points+hS/2);
    GI.PhiBgradS = GeneratePhiGrad1D(hS);
    GI.PhiBL = GeneratePhi1D(0,hL,hL,hL/2*points+hL/2);
    GI.PhiBgradL = GeneratePhiGrad1D(hL);
    
    %create element matrix
    GI.topology = 2;
    for i=1:ntot-1
        for j=1:GI.topology
            elmat(i,j) = i+j-1;
        end
    end
    
    %save all date to the structure
    GI.S = S;
    GI.L = L;
    GI.nS = nS;
    GI.nL = nL;
    GI.x = x;
    GI.ntot = ntot;
    GI.elmat = elmat;
    GI.hS = hS;
    GI.hL = hL;
end