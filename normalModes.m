% assume M = 1, Ze^2/(4pi epsilon_0) = 1
clear all
% global r r_d N omega_1_squared omega_z_squared;
[EigValues,EigVectors,wz,w1,positions,r_d,radius] = findNormalModes(4,3);
% save('NormalModes')
% save('NormalModes_nonspherical')
% clf
% hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
%              'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
%              'XLim',[-1.01 0.01],...               %#   set the x axis limit,
%              'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
%              'Color','none');               %#   and don't use a background color
% plot(EigValues,0,'bx');

% histogram(EigValues,50)

%find radii for N=6,8,30

% EigValueData = acrossAlphaData(3);

function [EigValueData,EigVectorData] = acrossAlphaData(N)
    alpha = [0.1 0.2 0.3 0.4 0.5 0.6 0.85 1 1.2 1.67 2 2.5 3.3 5 10];
    EigValueData = [];
    EigVectorData = [];
    for i = 1:length(alpha)
        [EigValues,EigVectors,wz,w1,r,r_d,radius]=findNormalModes(N,alpha(i));
        EigValueData = horzcat(EigValueData,EigValues);
        EigVectorData = cat(3,EigVectorData,EigVectors);
    end
end


function [EigValues,EigVectors,omega_z_squared,omega_1_squared,r,r_d,radius]=findNormalModes(N,alpha)
    global r r_d omega_z_squared omega_1_squared N;
    G = Greens(alpha);
    omega_z_squared = G;
    omega_1_squared = (1-G)/2;
%     omega_z_squared = 1;
%     omega_1_squared = 1;
    final_positions = crystal_graphs_energy(N,omega_z_squared,omega_1_squared);
    r = zeros(N,3);
    for i = 1:N
        r(i,:) = final_positions(i,1:3);
    end
    radius = sqrt(sum(r.^2,2));
    r_d = zeros(N,N);
    for i = 1:N
        for j = 1:N
            r_d(i,j) = sqrt(sum((r(i,:)-r(j,:)).^2));
        end
    end
    % r_d = (r_d + r_d') - eye(size(r_d,1)).*diag(r_d);

    D = zeros(3*N,3*N);

    for i = 1:N % x_i
        for j = 1:N % x_j
            % first N rows (x)
            D(i,j) = DXX(i,j);
            D(i,j+N) = DXY(i,j);
            D(i,j+2*N) = DXZ(i,j);
            % next N rows (y)
            D(i+N,j) = DYX(i,j);
            D(i+N,j+N) = DYY(i,j);
            D(i+N,j+2*N) = DYZ(i,j);
            % last N rows (z)
            D(i+2*N,j) = DZX(i,j);
            D(i+2*N,j+N) = DZY(i,j);
            D(i+2*N,j+2*N) = DZZ(i,j);
        end
    end
    [EigVectors, DiagValue] = eig(D);
    EigValues = diag(DiagValue);
end



function G = Greens(alpha)
    if alpha < 1
        beta = sqrt(1-alpha^2);
        G = 1/beta^2 - alpha/beta^3*atan(beta/alpha);
    elseif alpha > 1
        beta = sqrt(alpha^2-1);
        G = -1/beta^2 + alpha/(2*beta^3)*log(abs((alpha+beta)/(alpha-beta)));
    else
        G = 0.3334;
    end
end

function resXX = DXX(i,j)
    global r r_d N omega_1_squared;
    if i == j
        resXX = -omega_1_squared;
        for k = 1:N
            if i == k
                continue
            end
            resXX = resXX - (-1/r_d(i,k)^3 + 3/r_d(i,k)^5*(r(i,1)-r(k,1))^2);
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
            if i == k
                continue
            end
            resXY = resXY - 3/r_d(i,k)^5*(r(i,1)-r(k,1))*(r(i,2)-r(k,2));
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
            if i == k
                continue
            end
            resXZ = resXZ - 3/r_d(i,k)^5*(r(i,1)-r(k,1))*(r(i,3)-r(k,3));
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
            if i == k
                continue
            end
            resYY = resYY + (1/r_d(i,k)^3 - 3/r_d(i,k)^5*(r(i,2)-r(k,2))^2);
        end
    else
        resYY = - 1/r_d(i,j)^3 + 3/r_d(i,j)^5*(r(i,2)-r(j,2))^2;
    end
end
    
function resYX = DYX(i,j)
    global r r_d N;
    if i == j
        resYX = 0;
        for k = 1:N
            if i == k
                continue
            end
            resYX = resYX - 3/r_d(i,k)^5*(r(i,2)-r(k,2))*(r(i,1)-r(k,1));
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
            if i == k 
                continue, end
            resYZ = resYZ - 3/r_d(i,k)^5*(r(i,2)-r(k,2))*(r(i,3)-r(k,3));
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
            if i == k
                continue
            end
            resZZ = resZZ + 1/r_d(i,k)^3 - 3/r_d(i,k)^5*(r(i,3)-r(k,3))^2;
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
            if i == k
                continue
            end
            resZY = resZY - 3/r_d(i,k)^5*(r(i,3)-r(k,3))*(r(i,2)-r(k,2));
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
            if i == k
                continue
            end
            resZX = resZX - 3/r_d(i,k)^5*(r(i,3)-r(k,3))*(r(i,1)-r(k,1));
        end
    else
        resZX = 3/r_d(i,j)^5*(r(i,3)-r(j,3))*(r(i,1)-r(j,1));
    end
end