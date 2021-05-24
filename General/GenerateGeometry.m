%Generate the geometry information datastructure (GI) contain information about
%the finite element mesh and basis functions
function GI = GenerateGeometry(s,l,dn,ds,minvert,refineCenter)
    %make the empty datastructure
    GI = struct;
    
    %define the finite element mesh for the bridge
    poly1 = [3;4;-l/2;-s/2;-s/2;-l/2;ds;ds;0;0];
    poly2 = [3;4;s/2;l/2;l/2;s/2;ds;ds;0;0];
    poly3 = [3;4;-l/2;l/2;l/2;-l/2;0;0;-dn;-dn];
    gd = [poly1,poly2,poly3];
    ns = [80,80,80;49,50,51];
    sf = 'P1+P2+P3';
    Geometry = decsg(gd,sf,ns);
    [p,e,t] = initmesh(Geometry);

    %refine the mesh until the required number of vertices is achieved
    while length(p)<minvert
        [p,e,t] = refinemesh(Geometry,p,e,t);
    end
    
    %nodal positions in the different layers
    x = p(1,:); y = p(2,:);
    
    if refineCenter==true
        %needed for mesh refinement around corner
        it = find((mean(x(t(1:3,:)))>-1.1*s/2 & mean(x(t(1:3,:)))<1.1*s/2 & mean(y(t(1:3,:)))<dn/2 & mean(y(t(1:3,:)))>-dn));
        if isrow(it)
            it = it';
        end
        [p,e,t] = refinemesh(Geometry,p,e,t,it);
        x = p(1,:); y = p(2,:);
    end
    
    %position data in subdomains
    xs = x(y>=0);
    ys = y(y>=0);
    xn = x(y<=0);
    yn = y(y<=0);
    xs1 = xs(xs<=-0.9999*s/2);
    ys1 = ys(xs<=-0.9999*s/2);
    xs2 = xs(xs>=0.9999*s/2);
    ys2 = ys(xs>=0.9999*s/2);
    ns1 = length(xs1);
    ns2 = length(xs2);
    nn = length(xn);
    
    %element matrices of the different materials
    elmat = t(1:3,:);
    elmat = elmat';

    elmatn = elmat(mean(y(elmat),2)<=0,:);
    elmats = elmat(mean(y(elmat),2)>=0,:);
    elmats1 = elmats(mean(x(elmats),2)<=-s/2,:);
    elmats2 = elmats(mean(x(elmats),2)>=s/2,:);

    elmatbnd = e(1:2,:);
    elmatbnd = elmatbnd';

    elmatbnd1 = elmatbnd(sum(y(elmatbnd)==0,2)==2,:);
    elmatbnds1 = elmatbnd1(sum(x(elmatbnd1)<=-s/2,2)==2,:);
    elmatbnds2 = elmatbnd1(sum(x(elmatbnd1)>=s/2,2)==2,:);

    convs1 = zeros(length(x),1);
    convs1(unique(elmats1)) = 1:length(unique(elmats1));
    convs2 = zeros(length(x),1);
    convs2(unique(elmats2)) = 1:length(unique(elmats2));
    convn = zeros(length(x),1);
    convn(unique(elmatn)) = 1:length(unique(elmatn));
    
    %define element type by topology
    GI.topology = 3;
    GI.topologybnd = 2;
    
    %essential boundary condition points at the left and right edge of the superconductors
    counters1 = 0;
    counters2 = 0;
    boundarynodes = unique([e(1,:),e(2,:)]);
    for i = 1:length(boundarynodes)
        if x(boundarynodes(i)) == -l/2 && y(boundarynodes(i)) >=0
            counters1 = counters1+1;
            ess_bc_index_s1(counters1) = boundarynodes(i);
        end
        if x(boundarynodes(i)) == l/2 && y(boundarynodes(i)) >=0
            counters2 = counters2+1;
            ess_bc_index_s2(counters2) = boundarynodes(i);
        end     
    end
    GI.ess_bc_index_s1 = unique(ess_bc_index_s1);
    GI.ess_bc_index_s2 = unique(ess_bc_index_s2);
    
    GI.x = x;
    GI.y = y;
    GI.xn = xn;
    GI.yn = yn;
    GI.xs1 = xs1;
    GI.ys1 = ys1;
    GI.xs2 = xs2;
    GI.ys2 = ys2;
    GI.ns1 = ns1;
    GI.ns2 = ns2;
    GI.nn = nn;
    
    %element position matrices in different layers
    GI.elmats1 = elmats1;
    GI.elmats2 = elmats2;
    GI.elmatn = elmatn;
    GI.elmatbnds1 = elmatbnds1;
    GI.elmatbnds2 = elmatbnds2;
    
    %conversion between global elements and local elements
    GI.convs1 = convs1;
    GI.convs2 = convs2;
    GI.convn = convn; 
    GI.s = s;
    GI.l = l;
    GI.dn = dn;
    GI = CalculateElementData(GI); %calculate basis functions
end