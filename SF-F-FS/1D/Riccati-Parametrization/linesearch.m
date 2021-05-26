function alpha = linesearch(alpha_step,n_alpha,pk,f0,gammaF,gammaTF,GI,weights,omega,gammaBM,H,phase,Delta_0)
    
%initial choice of alpha;
alphalist = [];
min_alpha = 1;
flist = [];
%previous point
n = 0;
gammaF_old = gammaF;
gammaTF_old = gammaTF;
armijo = 10^(-4);
%compute at new point
while n<n_alpha
    alpha = min_alpha*(alpha_step^(n));
    alphalist = [alphalist,alpha];
    gammaF = gammaF_old-alphalist(end)*pk(1:GI.ntot);
    gammaTF = gammaTF_old-alphalist(end)*pk(GI.ntot+1:end);

    %determine new objective function value
    [S1,S1T,S2,S2T,S3,S3T,S4,S4T,S5,S5T,S6,S6T,f,fT] = BuildMatricesandVectors(GI,weights,omega,gammaBM,H,phase,gammaF,gammaTF,Delta_0);
    Stot = S1+S2+S3+S4+S5+S6;
    StotT = S1T+S2T+S3T+S4T+S5T+S6T;
    h = [Stot*gammaF-f;StotT*gammaTF-fT];
    
    flist = [flist,norm(h)];
    if flist(end)<(1-armijo*alpha)*f0
        return
    end
    n = n+1;
end
[~,indx] = min(flist);
alpha = alphalist(indx);
end