% assume M = 1, Ze^2/(4pi epsilon_0) = 1
clear all
global r r_d N;
alpha = 0.9;
N = 4;
G = Greens(alpha);
omega_z_squared = G;
omega_1_squared = (1-G)/2;
final_positions = crystal_graphs_energy(N,omega_z_squared,omega_1_squared);
r = zeros(4,3);
for i = 1:4
    r(i,:) = final_positions(i,1:3);
end
r_d = zeros(4,4);
for i = 1:4
    for j = 1:i
        r_d(i,j) = sqrt((r(i,1)-r(j,1))^2+(r(i,2)-r(j,2))^2+(r(i,3)-r(j,3))^2);
    end
end
r_d = (r_d + r_d') - eye(size(r_d,1)).*diag(r_d);

D = zeros(12,12);
% first N rows
for i = 1:N % x_i
    for j = 1:N % x_j
        if i == 1
            D(i,j) = -omega_1_squared + commonPartXX(i);
            D(i,j+4) = -commonPartXY(i);
            D(i,j+8) = -commonPartXZ(i);
        else
            D(i,j) = -commonPartXX(i);
            D(i,j+4) = commonPartXY(i);
            D(i,j+8) = commonPartXZ(i);
        end
    end
end

% next 4 rows
for i = 1:N % y_i
    for j = 1:N % y_j
        if i == 1
            D(N+i,j) = -commonPartYX(i);
            D(N+i,j+N) = -omega_1_squared + commonPartYY(i);
            D(N+i,j+2*N) = -commonPartYZ(i);
        else
            D(N+i,j) = commonPartYX(i);
            D(N+i,j+N) = -commonPartYY(i);
            D(N+i,j+2*N) = commonPartYZ(i);
        end
    end
end

% last four rows
for i = 1:N % z_i
    for j = 1:N %z_j
        if i == 1
            D(2*N+i,j) = -commonPartZX(i);
            D(2*N+i,j+N) = -commonPartZY(i);
            D(2*N+i,j+2*N) = -omega_z_squared + commonPartYZ(i);
        else
            D(2*N+i,j) = commonPartZX(i);
            D(2*N+i,j+N) = commonPartZY(i);
            D(2*N+i,j+2*N) = -commonPartZZ(i);
        end
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

function sumCommonPartXX = commonPartXX(ionIndex)
    global r r_d N;
    sumCommonPartXX = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartXX = sumCommonPartXX + 1/r_d(ionIndex,i)^3 - 3/r_d(ionIndex,i)^5*(r(ionIndex,1) - r(i,1))^2;
    end
end

function sumCommonPartXY = commonPartXY(ionIndex)
    global r r_d N;
    sumCommonPartXY = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartXY = sumCommonPartXY + 3/r_d(ionIndex,i)^5*(r(ionIndex,1) - r(i,1))*(r(ionIndex,2) - r(i,2));
    end
end

function sumCommonPartXZ = commonPartXZ(ionIndex)
    global r r_d N;
    sumCommonPartXZ = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartXZ = sumCommonPartXZ + 3/r_d(ionIndex,i)^5*(r(ionIndex,1) - r(i,1))*(r(ionIndex,3) - r(i,3));
    end
end

function sumCommonPartYY = commonPartYY(ionIndex)
    global r r_d N;
    sumCommonPartYY = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartYY = sumCommonPartYY + 1/r_d(ionIndex,i)^3 - 3/r_d(ionIndex,i)^5*(r(ionIndex,2) - r(i,2))^2;
    end
end

function sumCommonPartYX = commonPartYX(ionIndex)
    global r r_d N;
    sumCommonPartYX = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartYX = sumCommonPartYX + 3/r_d(ionIndex,i)^5*(r(ionIndex,1) - r(i,1))*(r(ionIndex,2) - r(i,2));
    end
end

function sumCommonPartYZ = commonPartYZ(ionIndex)
    global r r_d N;
    sumCommonPartYZ = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartYZ = sumCommonPartYZ + 3/r_d(ionIndex,i)^5*(r(ionIndex,2) - r(i,2))*(r(ionIndex,3) - r(i,3));
    end
end
    
function sumCommonPartZZ = commonPartZZ(ionIndex)
    global r r_d N;
    sumCommonPartZZ = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartZZ = sumCommonPartZZ + 1/r_d(ionIndex,i)^3 - 3/r_d(ionIndex,i)^5*(r(ionIndex,3) - r(i,3))^2;
    end
end

function sumCommonPartZY = commonPartZY(ionIndex)
    global r r_d N;
    sumCommonPartZY = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartZY = sumCommonPartZY + 3/r_d(ionIndex,i)^5*(r(ionIndex,3) - r(i,3))*(r(ionIndex,2) - r(i,2));
    end
end

function sumCommonPartZX = commonPartZX(ionIndex)
    global r r_d N;
    sumCommonPartZX = 0;
    for i = 1:N
        if i == ionIndex
            continue
        end
        sumCommonPartZX = sumCommonPartZX + 3/r_d(ionIndex,i)^5*(r(ionIndex,1) - r(i,1))*(r(ionIndex,3) - r(i,3));
    end
end