%Calculate the supercurrent through the N layer using the solution of the
%Usadel equations in the phi-parametrization
function Curr = CalculateCurrent(GI,n_mats,phi,T)

    %Approximates the gradient of the solution
    phigrad = zeros(GI.ntot,n_mats);
    for i=2:GI.nL-1
        phigrad(i,:) = (phi(i+1,:)-phi(i-1,:))/(2*GI.hL);
    end
    for i=GI.nL+1:GI.nL+GI.nS-2
        phigrad(i,:) = (phi(i+1,:)-phi(i-1,:))/(2*GI.hS);
    end
    for i=GI.nL+GI.nS:GI.ntot-1
        phigrad(i,:) = (phi(i+1,:)-phi(i-1,:))/(2*GI.hL);
    end

    %special treatment for gradient at center and at superconductor edge
    phigrad(1,:) = (phi(2,:)-phi(1,:))/(GI.hL);
    phigrad(GI.nL,:) = (phi(GI.nL+1,:)-phi(GI.nL,:))/(GI.hS);
    phigrad(GI.nL+GI.nS-1,:) = (phi(GI.nL+GI.nS-1,:)-phi(GI.nL+GI.nS-2,:))/(GI.hS);
    phigrad(GI.ntot,:) = (phi(GI.ntot,:)-phi(GI.ntot-1,:))/(GI.hL);
    
    %Calculate second Green's function
    G2 = zeros(GI.ntot,n_mats);
    for i=1:n_mats
        omega = (2*i-1)*T;
        G2(:,i) = 1./(omega^2+abs(phi(:,i)).^2);
    end
    
    %calculate the current
    Curr = (GI.S)*T*sum(imag(G2.*conj(phi).*phigrad),2);
end