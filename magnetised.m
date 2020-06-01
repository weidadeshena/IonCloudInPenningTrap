% assume we are in omega_c/2 frame, when the ion cloud is rotating in this
% frame, we say the ion cloud is magnetised.


N = 300;
omega_1_squared = 1;
omega_z_squared = 1;
positions = rand(N,3);
positions = update(N,positions,0,omega_1_squared,omega_z_squared,100);
positions = update(N,positions,0.1,omega_1_squared,omega_z_squared,100);
positions = update(N,positions,0.5,omega_1_squared,omega_z_squared,100);



function newPositions = update(N,oldPositions,omega_r_squared,omega_1_squared,omega_z_squared,steps)
% find the postions after a given number of steps
% omega_r squared is the angular velocity of the ion cloud in the rotating
% frame.
    stepSize = 0.01;
    theta = sqrt(omega_r_squared)*stepSize;
    R = [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
    s = scatter3(oldPositions(:,1),oldPositions(:,2),oldPositions(:,3),'x');
    axis([-10,10,-10,10,-10,10])
    for i = 1:steps
        F = force(N,oldPositions,omega_r_squared,omega_1_squared,omega_z_squared);
%         totalforce = sum(sum(F)) % print sumed forces
        newPositions = oldPositions + F*(0.5*(1/(1+exp(7.5*i/steps - 4)))); % proportional control
        % take rotation into account
        newPositions = (R*newPositions')';
%         % print rotation frequency
%         omega = angularFrequency(newPositions(:,1:2),oldPositions(:,1:2),stepSize)
        oldPositions = newPositions;
        set(s,'XData',newPositions(:,1),'YData',newPositions(:,2),'ZData',newPositions(:,3));
        drawnow
        pause(0.01)
    end
end

function omega = angularFrequency(newPositions,oldPositions,time)
% calculate the angualr velocity of the ion cloud
    dotProduct = sum(oldPositions.*newPositions,2);
    r1 = sqrt(sum(oldPositions.^2,2));
    r2 = sqrt(sum(newPositions.^2,2));
    theta = acos(dotProduct./(r1.*r2));
    omega = theta/time;
    omega = mean(omega);
end

% Cartesian Coord
% omega_r_squared is plasma rotation frequency in the rotating frame
% omega_1_squared is the radial trapping frequency
% omega_z_squared is the axial trapping frequency
function F = force(N,positions,omega_r_squared,omega_1_squared,omega_z_squared)
% calculate the force on each ion
    F = zeros(N,3);
    inverseDCubed = dist(positions,N);

    for i = 1:N
        % outwards
        % Coulomb
        for j = 1:N
            F(i,:) = F(i,:) + (positions(i,:) - positions(j,:))*inverseDCubed(i,j);
        end
        % centrifugal and trapping force
        F(i,1:2) = F(i,1:2) + (omega_r_squared - omega_1_squared) *positions(i,1:2);
        F(i,3) = F(i,3) - omega_z_squared * positions(i,3);
    end
end

function inverseDistanceCubed = dist(positions,N)
% calculate the distance between ion i and ion j, stored in a matrix
inverseDistance = inf(N,N);
for i = 1:N
    for j = i+1:N
        inverseDistance(i,j) = sqrt(sum((positions(i,:)-positions(j,:)).^2));
        inverseDistance(j,i) = inverseDistance(i,j);
    end
end
inverseDistanceCubed = inverseDistance.^-3;
end
    
