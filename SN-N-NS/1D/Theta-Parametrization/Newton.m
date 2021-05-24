%solves the nonlinear system from FEM using the Newton method
function [chi,theta,success] = Newton(GI,points,weights,phase,E,gamma,chi,theta,Delta_0,maxit,tol)
    success = true; %Does the Newton method converge?
    
    %iteration parameters
    Diff = 10;
    iter = 0;
    while ((Diff>tol) && (iter<maxit))
        [S1,f1,f2,f3,S2,f4,JS21,Jf11,Jf12,Jf21,Jf31,Jf32,Jf41,Jf42] = BuildMatricesandVectorsN(GI,weights,phase,E,gamma,chi,theta,Delta_0);
        
        %Build partial Jacobians and element vectors
        Jactheta1 = -Jf11-Jf21-Jf31;
        Jactheta2 = JS21;
        Jacchi1 = -Jf12-Jf32;
        Jacchi2 = sparse(GI.ntot,GI.ntot);
        f1tot = f1+f2+f3;
        f2tot = f4;

        %build full Jacobian and cost function
        Jactheta1 = Jactheta1+S1;
        Jactheta2 = Jactheta2-Jf41;
        Jacchi2 = Jacchi2+S2-Jf42;
        
        Jac = [Jactheta1,Jacchi1; Jactheta2,Jacchi2];
        h = [S1*theta-f1tot; S2*chi-f2tot];
                
        %find correct Newton direction and stepsize
        pk = Jac\h; %descent direction
        %apply linesearch
        alpha = linesearch(0.5,10,pk,norm(h),theta,chi,GI,weights,E,gamma,phase,Delta_0);
        
        theta = theta-alpha*pk(1:GI.ntot);
        chi = chi-alpha*pk(GI.ntot+1:end);
        
        Diff = norm(h);
        iter = iter+1;
        Difference(iter) = Diff;
    end
end