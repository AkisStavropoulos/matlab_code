

set1= s(1,1);

MasterPositionSet1= set1.master_position
MasterCartSet1= set1.master_cart

figure(1)
scatter(MasterCartSet1(:,3), MasterCartSet1(:,5))
axis equal