% attempts to find the structure of the inner crystal to justify why the
% density is the same. Failed.

% density_list = zeros(10);
% for i = 1:5
%     N = i*20+80;
%     density_list(i) = findCrystalDensity(N,3);
% end

% N = 500
% alpha = 0.3: ionCloud1
% alpha = 1: ionCloud2
% alpha = 3: ionCloud3
% N = 500;
% ionCloud1 = crystal_graphs_energy(N,0.6614,0.1693,0,0);
% ionCloud2 = crystal_graphs_energy(N,0.3334,0.3333,0,0);
% ionCloud3 = crystal_graphs_energy(N,0.1084,0.4456,0,0);
% ionCloud1 = ionCloud1(:,1:3);
% ionCloud2 = ionCloud2(:,1:3);
% ionCloud3 = ionCloud3(:,1:3);
% save 3ionCloud

innerSize = 3;
load 3ionCloud
innerCrystal1 = [];
innerCrystal2 = [];
innerCrystal3 = [];
for i = 1:N
    if all(abs(ionCloud1(i,:))<innerSize)
        innerCrystal1 = vertcat(innerCrystal1,ionCloud1(i,:));
    end
    if all(abs(ionCloud2(i,:))<innerSize)
        innerCrystal2 = vertcat(innerCrystal2,ionCloud2(i,:));
    end
    if all(abs(ionCloud3(i,:))<innerSize)
        innerCrystal3 = vertcat(innerCrystal3,ionCloud3(i,:));
    end
end

fig = figure;
subplot(1,3,1)
hold on
plot3(innerCrystal1(:,1),innerCrystal1(:,2),innerCrystal1(:,3),'b*','MarkerSize',4);

subplot(1,3,2)
plot3(innerCrystal2(:,1),innerCrystal2(:,2),innerCrystal2(:,3),'b*','MarkerSize',4);
subplot(1,3,3)
plot3(innerCrystal3(:,1),innerCrystal3(:,2),innerCrystal3(:,3),'b*','MarkerSize',4);

        

