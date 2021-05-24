function alpha = linesearch(alpha_step,n_alpha,pk,f0,GI,points,weights,points1D,weights1D,E,gamma_B,gamma,phase,ksi,Delta_0,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2)
%initial choice of alpha;
alphalist = [];
min_alpha = 1;
flist = [];
Difflist = [];
%previous point
n = 0;
armijo = 10^(-4);

thetas1_old = thetas1;
chis1_old = chis1;
thetas2_old = thetas2;
chis2_old = chis2;
thetan_old = thetan;
chin_old = chin;

%compute at new point
while n<n_alpha
    alpha = min_alpha*(alpha_step^(n));
    alphalist = [alphalist,alpha];

    thetas1 = thetas1_old-alphalist(end)*pk(1:GI.ns1);
    chis1 = chis1_old-alphalist(end)*pk(GI.ns1+1:2*GI.ns1);
    thetas2 = thetas2_old-alphalist(end)*pk(2*GI.ns1+1:2*GI.ns1+GI.ns2);
    chis2 = chis2_old-alphalist(end)*pk(2*GI.ns1+GI.ns2+1:2*(GI.ns1+GI.ns2));
    thetan = thetan_old-alphalist(end)*pk(2*(GI.ns1+GI.ns2)+1:2*(GI.ns1+GI.ns2)+GI.nn);
    chin = chin_old-alphalist(end)*pk(2*(GI.ns1+GI.ns2)+GI.nn+1:end);

    %determine new objective function value
    [S11,S21,f11,f21,f31,f41,fB11,fB21,S12,S22,f12,f22,f32,f42,fB12,fB22,S1n,S2n,f1n,f2n,fB1n,fB2n] = BuildVecMatricesandVectors(GI,points,weights,points1D,weights1D,gamma_B,gamma,E,ksi,thetas1,chis1,thetas2,chis2,thetan,chin,deltas1,deltas2);
    Stot = blkdiag(S11,S21,S12,S22,S1n,S2n);
    ftot = [f11+f21+f31+fB11;f41+fB21;f12+f22+f32+fB12;f42+fB22;f1n+f2n+fB1n;fB2n];
    [~,Stot,ftot] = ProcessEssentialBoundaryConditions(sparse(2*(GI.ns1+GI.ns2+GI.nn),2*(GI.ns1+GI.ns2+GI.nn)),Stot,ftot,GI,phase,Delta_0,E);
    h = Stot*[thetas1;chis1;thetas2;chis2;thetan;chin]-ftot;
    
    flist = [flist,norm(h)];
    if flist(end)<(1-armijo*alpha)*f0
        return
    end
    n = n+1;
end
end