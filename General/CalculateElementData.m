%Calculate linear basis function data for all triangular elements in FEM grid
function GI = CalculateElementData(GI)
    
    %geometry data of the first superconductor
    GI.Deltas1 = zeros(length(GI.elmats1(:,1)),1);
    GI.as1 = zeros(length(GI.elmats1(:,1)),3);
    GI.bs1 = zeros(length(GI.elmats1(:,1)),3);
    GI.cs1 = zeros(length(GI.elmats1(:,1)),3);
    %geometry data of the second superconductor
    GI.Deltas2 = zeros(length(GI.elmats2(:,1)),1);
    GI.as2 = zeros(length(GI.elmats2(:,1)),3);
    GI.bs2 = zeros(length(GI.elmats2(:,1)),3);
    GI.cs2 = zeros(length(GI.elmats2(:,1)),3);
    %geometry data of the normal metal
    GI.Deltan = zeros(length(GI.elmatn(:,1)),1);
    GI.an = zeros(length(GI.elmatn(:,1)),3);
    GI.bn = zeros(length(GI.elmatn(:,1)),3);
    GI.cn = zeros(length(GI.elmatn(:,1)),3);
    
    %caculate element data in first superconductor
    for i=1:length(GI.elmats1)
        for index1 = 1:GI.topology
            xc(index1) = GI.x(GI.elmats1(i,index1));
            yc(index1) = GI.y(GI.elmats1(i,index1));
        end
        GI.Deltas1(i) = xc(2)*yc(3)-xc(3)*yc(2)+xc(3)*yc(1)-xc(1)*yc(3)+xc(1)*yc(2)-xc(2)*yc(1);
        B_mat = [1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)] \  eye(3);

        GI.as1(i,:) = B_mat(1,1:3);
        GI.bs1(i,:) = B_mat(2,1:3);
        GI.cs1(i,:) = B_mat(3,1:3);
    end
    %caculate element data in second superconductor
    for i=1:length(GI.elmats2)
        for index1 = 1:GI.topology
            xc(index1) = GI.x(GI.elmats2(i,index1));
            yc(index1) = GI.y(GI.elmats2(i,index1));
        end
        GI.Deltas2(i) = xc(2)*yc(3)-xc(3)*yc(2)+xc(3)*yc(1)-xc(1)*yc(3)+xc(1)*yc(2)-xc(2)*yc(1);
        B_mat = [1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)] \  eye(3);

        GI.as2(i,:) = B_mat(1,1:3);
        GI.bs2(i,:) = B_mat(2,1:3);
        GI.cs2(i,:) = B_mat(3,1:3);
    end
    %caculate element data in normal metal
    for i=1:length(GI.elmatn)
        for index1 = 1:GI.topology
            xc(index1) = GI.x(GI.elmatn(i,index1));
            yc(index1) = GI.y(GI.elmatn(i,index1));
        end
        GI.Deltan(i) = xc(2)*yc(3)-xc(3)*yc(2)+xc(3)*yc(1)-xc(1)*yc(3)+xc(1)*yc(2)-xc(2)*yc(1);
        B_mat = [1 xc(1) yc(1);1 xc(2) yc(2);1 xc(3) yc(3)] \  eye(3);

        GI.an(i,:) = B_mat(1,1:3);
        GI.bn(i,:) = B_mat(2,1:3);
        GI.cn(i,:) = B_mat(3,1:3);
    end
end