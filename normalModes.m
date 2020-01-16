% assume M = 1, Ze^2/(4pi epsilon_0) = 1
clear all
global r r_d N omega_1_squared omega_z_squared;
alpha = 0.4;
N = 10;
G = Greens(alpha);
omega_z_squared = G;
omega_1_squared = (1-G)/2;
final_positions = crystal_graphs_energy(N,omega_z_squared,omega_1_squared);
r = final_positions(:,1:3);
r_d = zeros(N,N);
for i = 1:N
    for j = 1:N
        r_d(i,j) = sqrt(sum((r(i,:)-r(j,:)).^2));
    end
end
r_d = (r_d + r_d') - eye(size(r_d,1)).*diag(r_d);

D = zeros(3*N,3*N);
% first N rows
for i = 1:N % x_i
    for j = 1:N % x_j
        D(i,j) = DXX(i,j);
        D(i,j+N) = DXY(i,j);
        D(i,j+2*N) = DXZ(i,j);
    end
end
% next 4 rows
for i = 1:N % y_i
    for j = 1:N % y_j
        D(i+N,j) = DYX(i,j);
        D(i+N,j+N) = DYY(i,j);
        D(i+N,j+2*N) = DYZ(i,j);
    end
end
% last 4 rows
for i = 1:N % z_i
    for j = 1:N %z_j
        D(i+2*N,j) = DZX(i,j);
        D(i+2*N,j+N) = DZY(i,j);
        D(i+2*N,j+2*N) = DZZ(i,j);
    end
end

[EigVectors, DiagValue] = eig(D);
EigValues = diag(DiagValue);

function G = Greens(alpha)
    if alpha < 1
        beta = sqrt(1-alpha^2);
        G = 1/beta^2 - alpha/beta^3*atan(beta/alpha);
    elseif alpha > 1
        beta = sqrt(alpha^2-1);
        G = -1/beta^2 + alpha/(2*beta^3)*log(abs((alpha+beta)/(alpha-beta)));
    end
end

function resXX = DXX(i,j)
    global r r_d N omega_1_squared;
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

function resXY = DXY(i,j)
    global r r_d N;
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

function resXZ = DXZ(i,j)
    global r r_d N;
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

function resYY = DYY(i,j)
    global r r_d N omega_1_squared;
    if i == j
        resYY = -omega_1_squared;
        for k = 1:N
            if i ~= k
                resYY = resYY + (1/r_d(i,k)^3 - 3/r_d(i,k)^5*(r(i,2)-r(k,2))^2);
            end
         end
    else
        resYY = - (1/r_d(i,j)^3 + 3/r_d(i,j)^5*(r(i,2)-r(j,2))^2);
    end
end
    
function resYX = DYX(i,j)
    global r r_d N;
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

function resYZ = DYZ(i,j)
    global r r_d N;
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
    
function resZZ = DZZ(i,j)
    global r r_d N omega_z_squared;
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

function resZY = DZY(i,j)
    global r r_d N;
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

function resZX = DZX(i,j)
    global r r_d N;
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