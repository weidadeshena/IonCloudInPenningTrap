% find the density of an ion cloud
N = 100;
alpha = 0.3;
findDensity(N,alpha)

function density = findDensity(N,alpha)
% find the steady configuration of the ion cloud first, then circle some
% area inside the ion cloud, count the number of ions in the selected
% sphere, then divided by the volume of sphere. The ion cloud is sampled
% 125 times
[om_z_squared,om_1_squared] = Greens(alpha);
positions = crystal_graphs_energy(N,om_z_squared,om_1_squared,0,0);
% find the "radius" for the x,y,z dimension
ionRange = floor((max(positions) - min(positions))/2);
% set the sample radius to be 1/2 of the smallest dimenstion
sphereRadius = min(ionRange)*0.5;
% set the increment to be proportion to the "radius" of that dimension
xInc = ionRange(1)*0.2;
yInc = ionRange(2)*0.2;
zInc = ionRange(3)*0.2;

% start counting the number of ions in the sampled sphere
ionCounter = 0;
for i = -2:2
    for j = -2:2
        for k = -2:2
            ionCounter = ionCounter + inSphere(i*xInc,j*yInc,k*zInc,sphereRadius,positions);
        end
    end
end
% find the averge number per sample sphere
numberPerSphere = ionCounter/(5^3);
% find the density
density = numberPerSphere/(4*pi/3*sphereRadius^3);

end

function numberOfIons = inSphere(sphereX,sphereY,sphereZ,sphereR,positions)
% find the number of ions in a sphere
    numberOfIons = 0;
    % move the frame so that the ion is in origin and the cube is in ion
    % frame
    for i = 1:length(positions)
        ionInSphereFrame = [positions(i,1)-sphereX positions(i,2)-sphereY positions(i,3)-sphereZ];
        if sum(ionInSphereFrame.^2) < sphereR^2
            numberOfIons = numberOfIons + 1;
        end
    end
end
