function Curr = CalculateCurrent(GI,n_mats,phiP,phiCM,T,H)
    %calculate the gradient of the solutions
    phigrad_p = zeros(GI.ntot,n_mats);
    phigrad_m = zeros(GI.ntot,n_mats);
    for i=2:GI.nL-1
        phigrad_p(i,:) = (phiP(i+1,:)-phiP(i-1,:))/(2*GI.hL);
        phigrad_m(i,:) = (phiCM(i+1,:)-phiCM(i-1,:))/(2*GI.hL);
    end
    for i=GI.nL+1:GI.nL+GI.nS-2
        phigrad_p(i,:) = (phiP(i+1,:)-phiP(i-1,:))/(2*GI.hS);
        phigrad_m(i,:) = (phiCM(i+1,:)-phiCM(i-1,:))/(2*GI.hS);
    end
    for i=GI.nL+GI.nS:GI.ntot-1
        phigrad_p(i,:) = (phiP(i+1,:)-phiP(i-1,:))/(2*GI.hL);
        phigrad_m(i,:) = (phiCM(i+1,:)-phiCM(i-1,:))/(2*GI.hL);
    end
    %special treatment for gradient at superconductor edge
    phigrad_p(1,:) = (phiP(2,:)-phiP(1,:))/(GI.hL);
    phigrad_m(1,:) = (phiCM(2,:)-phiCM(1,:))/(GI.hL);
    phigrad_p(GI.nL,:) = (phiP(GI.nL+1,:)-phiP(GI.nL,:))/(GI.hS);
    phigrad_m(GI.nL,:) = (phiCM(GI.nL+1,:)-phiCM(GI.nL,:))/(GI.hS);
    phigrad_p(GI.nL+GI.nS-1,:) = (phiP(GI.nL+GI.nS-1,:)-phiP(GI.nL+GI.nS-2,:))/(GI.hS);
    phigrad_m(GI.nL+GI.nS-1,:) = (phiCM(GI.nL+GI.nS-1,:)-phiCM(GI.nL+GI.nS-2,:))/(GI.hS);
    phigrad_p(GI.ntot,:) = (phiP(GI.ntot,:)-phiP(GI.ntot-1,:))/(GI.hL);
    phigrad_m(GI.ntot,:) = (phiCM(GI.ntot,:)-phiCM(GI.ntot-1,:))/(GI.hL);
    
    %Calculate regular Green's function
    G2_p = zeros(GI.ntot,n_mats);
    G2_m = zeros(GI.ntot,n_mats);
    for i=1:n_mats
        omega = (2*i-1)*T;
        omega_p = omega+1i*H;
        omega_m = -omega+1i*H;
        G2_p(:,i) = 1./(omega_p^2+phiP(:,i).*phiCM(:,i));
        G2_m(:,i) = 1./(omega_m^2+conj(phiCM(:,i)).*conj(phiP(:,i)));
    end
    
    %calculate the current
    Curr = 1i/4*(GI.S)*T*sum(G2_p.*(phiP.*phigrad_m-phiCM.*phigrad_p),2);
    Curr = Curr+1i/4*(GI.S)*T*sum(G2_m.*(conj(phiCM).*conj(phigrad_p)-conj(phiP).*conj(phigrad_m)),2);
end