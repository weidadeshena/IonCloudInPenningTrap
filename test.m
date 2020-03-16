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