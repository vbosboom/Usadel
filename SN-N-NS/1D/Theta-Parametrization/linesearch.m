%Performs line search on the descent direction pk to find the optimal step
%size alpha
function alpha = linesearch(alpha_step,n_alpha,pk,f0,theta,chi,GI,weights,E,gamma,phase,Delta_0)

%initial choice of alpha;
alphalist = [];
min_alpha = 1;
flist = [];
%previous point
n = 0;
theta_old = theta;
chi_old = chi;
armijo = 10^(-4);
%compute at new point
while n<n_alpha
    alpha = min_alpha*(alpha_step^(n));
    alphalist = [alphalist,alpha];
    theta = theta_old-alphalist(end)*pk(1:GI.ntot);
    chi = chi_old-alphalist(end)*pk(GI.ntot+1:end);

    %determine new objective function value
    [S1,f1,f2,f3,S2,f4] = BuildMatricesandVectors(GI,weights,phase,E,gamma,chi,theta,Delta_0);
    f1tot = f1+f2+f3;
    f2tot = f4;
    h = [S1*theta-f1tot; S2*chi-f2tot];
    
    flist = [flist,norm(h)];
    if flist(end)<(1-armijo*alpha)*f0
        return
    end
    n = n+1;
end
[~,indx] = min(flist);
alpha = alphalist(indx);
end