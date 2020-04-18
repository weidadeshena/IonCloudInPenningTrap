
% density_list = zeros(10);
% for i = 1:5
%     N = i*20+80;
%     density_list(i) = findCrystalDensity(N,3);
% end

% N = 500
% alpha = 0.3: ionCloud1
% alpha = 1: ionCloud2
% alpha = 3: ionCloud3
N = 500;
ionCloud1 = crystal_graphs_energy(N,0.6614,0.1693,0,0);
ionCloud2 = crystal_graphs_energy(N,0.3334,0.3333,0,0);
ionCloud3 = crystal_graphs_energy(N,0.1084,0.4456,0,0);

for i = 1:N
    if any(ionCloud1(i,:)>3)
        ionCloud1(i,:) = [];
    end
    if any(ionCloud2(i,:)>3)
        ionCloud2(i,:) = [];
    end
    if any(ionCloud3(i,:)>3)
        ionCloud3(i,:) = [];
    end
end

fig = figure;
subplot(1,3,1)
plot3(ionCloud1(:,1),ionCloud1(:,2),ionCloud1(:,3),'b*','MarkerSize',4);
subplot(1,3,2)
plot3(ionCloud2(:,1),ionCloud2(:,2),ionCloud2(:,3),'b*','MarkerSize',4);
subplot(1,3,3)
plot3(ionCloud3(:,1),ionCloud3(:,2),ionCloud3(:,3),'b*','MarkerSize',4);

        

k = 20;
N = 50;
separation = zeros(k,1);
omega_z = 0.334;
omega_1 = 0.333;
omega_r = linspace(-omega_1+0.01,omega_1-0.01,k);
for i = 1:k
    positions = crystal_graphs_energy(N,omega_z,omega_1-omega_r(i),0,0);
    distance = dist(positions,N);
    separation(i) = mean(min(distance));
end


function Distance = dist(positions,N)
Distance = inf(N,N);
for i = 1:N
    for j = i+1:N
        Distance(i,j) = sqrt(sum((positions(i,:)-positions(j,:)).^2));
        Distance(j,i) = Distance(i,j);
    end
end
end

