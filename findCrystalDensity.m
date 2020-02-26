function numberDensity = findCrystalDensity(N,alpha)
[om_z_squared,om_1_squared] = Greens(alpha);
positions = crystal_graphs_energy(N,om_z_squared,om_1_squared,0,0);

% load NormalModes_1000 positions
ionRange = floor((max(positions) - min(positions))/2);
sphereRadius = min(ionRange)*0.5;
xInc = ionRange(1)*0.2;
yInc = ionRange(2)*0.2;
zInc = ionRange(3)*0.2;

ionCounter = 0;
for i = -2:2
    for j = -2:2
        for k = -2:2
            ionCounter = ionCounter + inSphere(i*xInc,j*yInc,k*zInc,sphereRadius,positions);
        end
    end
end

numberPerSphere = ionCounter/(5^3);
numberDensity = numberPerSphere/(4*pi/3*sphereRadius^3);
    
end

function numberOfIons = inSphere(sphereX,sphereY,sphereZ,sphereR,positions)
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
