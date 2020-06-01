load NormalModes positions r_d wz
N = 4;
% my attempts for guiding centre approximation.... didn't work
timestep = 1000;
dt = 0.01;
B = 2;
U_old = 0;
positions = positions';

for i = 1:timestep
    [newPositions,U] = dxbydt(positions,r_d,N,wz,dt,U_old,B);
    positions = newPositions;
    U_old = U;
end
positions
U

function [newPositions,U] = dxbydt(positions,r_d,N,omega_z_squared,dt,U_old,B)
[E,U] = dUbydt(positions,r_d,N,omega_z_squared,dt,U_old);
newPositions = zeros(3,N);
for i = 1:N
    newPositions(:,i) = (cross(E(:,:,i),B*[0;0;1])/B^2 + U(:,:,i))*dt;
end
newPositions = newPositions + positions;
end

function [E,U] = dUbydt(positions,r_d,N,omega_z_squared,dt,U_old)
e = 1;
U = zeros(3,1,N);
E = E_field(positions,r_d,N,omega_z_squared);
U = e*E*dt;
U = U + U_old;
end

function E = E_field(positions,r_d,N,omega_z_squared)
E = zeros(3,1,N);
for i = 1:N
    for j = 1:N
        if i~=j
            E(:,1,i) = E(:,1,i) + r_d(i,j)^-1.5*(positions(:,i) - positions(:,j));
        end
    end
    E(:,1,i) = E(:,1,i) + [0.5 0.5 -1]*omega_z_squared*positions(:,i);
end
end