%Calculates the BCS energy gap at temperature T normalized 
%to pi*k_b*T_c according to http://ltl.tkk.fi/research/theory/qc/bcsgap.pdf
function Delta = BCSGap(T)
m = 1000;
Delta = 0.1;
Diff = 10;
tol = 10^(-8);
while Diff>tol
    Delta_old = Delta;
    Sumf = 0;
    SumOmega=0;
    SumfDer = 0;
    for i=1:m
        omega = i-1/2;
        SumOmega = SumOmega+1/omega;
        Sumf = Sumf+T/sqrt(T^2*omega^2+Delta^2);
        SumfDer = SumfDer+T*Delta/((T^2*omega^2+Delta^2)^(3/2));
    end
    G = SumOmega-Sumf+log((m*T+sqrt(m^2*T^2+Delta^2))/(2*m))-1/24*(1/m^2-m*T^3/(m^2*T^2+Delta^2)^(3/2));
    Gdiff = SumfDer+Delta/(sqrt(m^2*T^2+Delta^2)*(m*T+sqrt(m^2*T^2+Delta^2))) - m*T^3*Delta/(8*(m^2*T^2+Delta^2))^(5/2);
    Delta = Delta-G/Gdiff;
    Diff = norm(Delta-Delta_old);
end
Delta = Delta*2;
end