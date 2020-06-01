function [EigValues,EigVectors,r]=findNormalModes(N,alpha)
% analytically find the normal modes. the eigenvalues are -frequency^2 and
% the eigen vectors are the direction of the eigenmodes
    [omega_z_squared,omega_1_squared] = Greens(alpha);
    final_positions = crystal_graphs_energy(N,omega_z_squared,omega_1_squared,0,0);
    r = zeros(N,3);
    for i = 1:N
        r(i,:) = final_positions(i,1:3);
    end
    r_d = zeros(N,N);
    for i = 1:N
        for j = 1:N
            r_d(i,j) = sqrt(sum((r(i,:)-r(j,:)).^2));
        end
    end

    D = zeros(3*N,3*N);

    for i = 1:N % x_i
        for j = 1:N % x_j
            % first N rows (x)
            D(i,j) = DXX(i,j,N,r,r_d,omega_1_squared);
            D(i,j+N) = DXY(i,j,N,r,r_d);
            D(i,j+2*N) = DXZ(i,j,N,r,r_d);
            % next N rows (y)
            D(i+N,j) = DYX(i,j,N,r,r_d);
            D(i+N,j+N) = DYY(i,j,N,r,r_d,omega_1_squared);
            D(i+N,j+2*N) = DYZ(i,j,N,r,r_d);
            % last N rows (z)
            D(i+2*N,j) = DZX(i,j,N,r,r_d);
            D(i+2*N,j+N) = DZY(i,j,N,r,r_d);
            D(i+2*N,j+2*N) = DZZ(i,j,N,r,r_d,omega_z_squared);
        end
    end
    [EigVectors, DiagValue] = eig(D);
    EigValues = diag(DiagValue);
end

function resXX = DXX(i,j,N,r,r_d,omega_1_squared)
    if i == j
        resXX = -omega_1_squared;
        for k = 1:N
            if i ~= k
                resXX = resXX - (-1/r_d(i,k)^3 + 3/r_d(i,k)^5*(r(i,1)-r(k,1))^2);
            end
        end
    else
        resXX = - 1/r_d(i,j)^3 + 3/r_d(i,j)^5*(r(i,1) - r(j,1))^2;
    end
end

function resXY = DXY(i,j,N,r,r_d)
    if i == j
        resXY = 0;
        for k = 1:N
            if i ~= k
                resXY = resXY - 3/r_d(i,k)^5*(r(i,1)-r(k,1))*(r(i,2)-r(k,2));
            end
        end
    else       
        resXY = 3/r_d(i,j)^5*(r(i,1)-r(j,1))*(r(i,2)-r(j,2));
    end
end

function resXZ = DXZ(i,j,N,r,r_d)
    if i==j
        resXZ = 0;
        for k = 1:N
            if i ~= k
                resXZ = resXZ - 3/r_d(i,k)^5*(r(i,1)-r(k,1))*(r(i,3)-r(k,3));
            end
         end
    else
        resXZ = 3/r_d(i,j)^5*(r(i,1)-r(j,1))*(r(i,3)-r(j,3));
    end
end

function resYY = DYY(i,j,N,r,r_d,omega_1_squared)
    if i == j
        resYY = -omega_1_squared;
        for k = 1:N
            if i ~= k
                resYY = resYY + (1/r_d(i,k)^3 - 3/r_d(i,k)^5*(r(i,2)-r(k,2))^2);
            end
         end
    else
        resYY = - 1/r_d(i,j)^3 + 3/r_d(i,j)^5*(r(i,2)-r(j,2))^2;
    end
end
    
function resYX = DYX(i,j,N,r,r_d)
    if i == j
        resYX = 0;
        for k = 1:N
            if i ~= k
                resYX = resYX - 3/r_d(i,k)^5*(r(i,2)-r(k,2))*(r(i,1)-r(k,1));
            end
        end
    else
        resYX = 3/r_d(i,j)^5*(r(i,2)-r(j,2))*(r(i,1)-r(j,1));
    end
end

function resYZ = DYZ(i,j,N,r,r_d)
    if i==j
        resYZ = 0;
        for k = 1:N
            if i ~= k 
                resYZ = resYZ - 3/r_d(i,k)^5*(r(i,2)-r(k,2))*(r(i,3)-r(k,3));
            end
        end
    else
        resYZ = 3/r_d(i,j)^5*(r(i,2)-r(j,2))*(r(i,3)-r(j,3));
    end
end
    
function resZZ = DZZ(i,j,N,r,r_d,omega_z_squared)
    if i == j
        resZZ = -omega_z_squared;
        for k = 1:N
            if i ~= k
                resZZ = resZZ + 1/r_d(i,k)^3 - 3/r_d(i,k)^5*(r(i,3)-r(k,3))^2;
            end
        end
    else
        resZZ = -1/r_d(i,j)^3 + 3/r_d(i,j)^5*(r(i,3)-r(j,3))^2;
    end
end

function resZY = DZY(i,j,N,r,r_d)
    if i == j
        resZY = 0;
        for k = 1:N
            if i ~= k
                resZY = resZY - 3/r_d(i,k)^5*(r(i,3)-r(k,3))*(r(i,2)-r(k,2));
            end
        end
    else
        resZY = 3/r_d(i,j)^5*(r(i,3)-r(j,3))*(r(i,2)-r(j,2));
    end
end

function resZX = DZX(i,j,N,r,r_d)
    if i == j
        resZX = 0;
        for k = 1:N
            if i ~= k
                resZX = resZX - 3/r_d(i,k)^5*(r(i,3)-r(k,3))*(r(i,1)-r(k,1));
            end
        end
    else
        resZX = 3/r_d(i,j)^5*(r(i,3)-r(j,3))*(r(i,1)-r(j,1));
    end
end