function [Currxs1,Currxs2,Currxn,Currys1,Currys2,Curryn,indx] = CalculateCurrent(GI,n_mats,phis1,phis2,phin,T,gamma_N,gamma_S,ksi)
%     hold on
    %for the first superconductor
    coefgradxs1 = zeros(length(GI.elmats1),n_mats);
    coefgradys1 = zeros(length(GI.elmats1),n_mats);
    coefs1 = zeros(length(GI.elmats1),n_mats);
    G2s1 = zeros(length(GI.elmats1),n_mats);
    xmids1 = zeros(length(GI.elmats1),1);
    ymids1 = zeros(length(GI.elmats1),1);
    
    %calculate functions and gradients in center of elements
    for i=1:length(GI.elmats1)
        xmids1(i) = mean(GI.x(GI.elmats1(i,:)));
        ymids1(i) = mean(GI.y(GI.elmats1(i,:)));
        for index1 = 1:3
            coefs1(i,:) = coefs1(i,:)+phis1(GI.convs1(GI.elmats1(i,index1)),:)*(GI.as1(i,index1)+GI.bs1(i,index1)*xmids1(i)+GI.cs1(i,index1)*ymids1(i));
            coefgradxs1(i,:) = coefgradxs1(i,:)+phis1(GI.convs1(GI.elmats1(i,index1)),:)*GI.bs1(i,index1);
            coefgradys1(i,:) = coefgradys1(i,:)+phis1(GI.convs1(GI.elmats1(i,index1)),:)*GI.cs1(i,index1);
        end
    end
    
    for i = 1:n_mats
        omega = (2*i-1)*T;
        G2s1(:,i) = 1./(omega^2+abs(coefs1(:,i)).^2);
    end 
    CurrFacResist = GI.s*(1+2*sqrt(gamma_N)/(GI.s)*coth((GI.l-GI.s)/(2*sqrt(gamma_N))));
    CurrFacPairBreak = (1/0.385)*ksi*(21*sqrt(3)*zeta(3)/(4*pi));
    Currxs1 = T*sum(imag(G2s1.*conj(coefs1).*coefgradxs1),2)*CurrFacPairBreak;
    Currys1 = T*sum(imag(G2s1.*conj(coefs1).*coefgradys1),2)*CurrFacPairBreak;
%     fill3(GI.x(GI.elmats1)',GI.y(GI.elmats1)',kron(ones(1,size(GI.elmats1,2)),Currxs1)',kron(ones(1,size(GI.elmats1,2)),Currxs1)');
    
    %for the second superconductor
    coefgradxs2 = zeros(length(GI.elmats2),n_mats);
    coefgradys2 = zeros(length(GI.elmats2),n_mats);
    coefs2 = zeros(length(GI.elmats2),n_mats);
    G2s2 = zeros(length(GI.elmats2),n_mats);
    xmids2 = zeros(length(GI.elmats2),1);
    ymids2 = zeros(length(GI.elmats2),1);
    for i=1:length(GI.elmats2)
        xmids2(i) = mean(GI.x(GI.elmats2(i,:)));
        ymids2(i) = mean(GI.y(GI.elmats2(i,:)));
        for index1 = 1:3
            coefs2(i,:) = coefs2(i,:)+phis2(GI.convs2(GI.elmats2(i,index1)),:)*(GI.as2(i,index1)+GI.bs2(i,index1)*xmids2(i)+GI.cs2(i,index1)*ymids2(i));
            coefgradxs2(i,:) = coefgradxs2(i,:)+phis2(GI.convs2(GI.elmats2(i,index1)),:)*GI.bs2(i,index1);
            coefgradys2(i,:) = coefgradys2(i,:)+phis2(GI.convs2(GI.elmats2(i,index1)),:)*GI.cs2(i,index1);
        end
    end
    for i = 1:n_mats
        omega = (2*i-1)*T;
        G2s2(:,i) = 1./(omega^2+abs(coefs2(:,i)).^2);
    end
    CurrFacResist = GI.s*(1+2*sqrt(gamma_N)/(GI.s)*coth((GI.l-GI.s)/(2*sqrt(gamma_N))));
    CurrFacPairBreak = (1/0.385)*ksi*(21*sqrt(3)*zeta(3)/(4*pi));
    Currxs2 = T*sum(imag(G2s2.*conj(coefs2).*coefgradxs2),2)*CurrFacPairBreak;
    Currys2 = T*sum(imag(G2s2.*conj(coefs2).*coefgradys2),2)*CurrFacPairBreak;
%     fill3(GI.x(GI.elmats2)',GI.y(GI.elmats2)',kron(ones(1,size(GI.elmats2,2)),Currxs2)',kron(ones(1,size(GI.elmats2,2)),Currxs2)');
    
    %for the normal metal
    coefgradxn = zeros(length(GI.elmatn),n_mats);
    coefgradyn = zeros(length(GI.elmatn),n_mats);
    coefn = zeros(length(GI.elmatn),n_mats);
    G2n = zeros(length(GI.elmatn),n_mats);
    xmidn = zeros(length(GI.elmatn),1);
    ymidn = zeros(length(GI.elmatn),1);
    
    for i=1:length(GI.elmatn)
        xmidn(i) = mean(GI.x(GI.elmatn(i,:)));
        ymidn(i) = mean(GI.y(GI.elmatn(i,:)));
        for index1 = 1:3
            coefn(i,:) = coefn(i,:)+phin(GI.convn(GI.elmatn(i,index1)),:)*(GI.an(i,index1)+GI.bn(i,index1)*xmidn(i)+GI.cn(i,index1)*ymidn(i));
            coefgradxn(i,:) = coefgradxn(i,:)+phin(GI.convn(GI.elmatn(i,index1)),:)*GI.bn(i,index1);
            coefgradyn(i,:) = coefgradyn(i,:)+phin(GI.convn(GI.elmatn(i,index1)),:)*GI.cn(i,index1);
        end
    end
    
    for i = 1:n_mats
        omega = (2*i-1)*T;
        G2n(:,i) = 1./(omega^2+abs(coefn(:,i)).^2);
    end
    [~,indx] = min(abs(xmidn-0).^2+abs(ymidn+GI.dn/2).^2);
    CurrFacResist = gamma_N/gamma_S/ksi*GI.s*(1+2*sqrt(gamma_N)/(GI.s)*coth((GI.l-GI.s)/(2*sqrt(gamma_N))));
    CurrFacPairBreak = (1/0.385)*gamma_N/gamma_S*(21*sqrt(3)*zeta(3)/(4*pi));
%     Currxn = GI.s*T*sum(imag(G2n.conj(coefn).*coefgradxn),2)(1+2*sqrt(gamma_N)/(GI.s)*coth((GI.l-GI.s)/(2*sqrt(gamma_N))));
%     Curryn = GI.s*T*sum(imag(G2n.conj(coefn).*coefgradyn),2)(1+2*sqrt(gamma_N)/(GI.s)*coth((GI.l-GI.s)/(2*sqrt(gamma_N))));
    Currxn = T*sum(imag(G2n.*conj(coefn).*coefgradxn),2)*CurrFacPairBreak;
    Curryn = T*sum(imag(G2n.*conj(coefn).*coefgradyn),2)*CurrFacPairBreak;%     figure
%     fill3(GI.x(GI.elmatn)',GI.y(GI.elmatn)',kron(ones(1,size(GI.elmatn,2)),Currxn)',kron(ones(1,size(GI.elmatn,2)),Currxn)');
%     colorbar
%     title('current density in the x-direction')
end