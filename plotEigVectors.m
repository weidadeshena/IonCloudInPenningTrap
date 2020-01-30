% clear all

plotMode(6)

function plotMode(j)
    load NormalModes EigValues EigVectors positions
% global EigValues EigVectors positions;
    axisMin = min(positions(:))-0.1;
    axisMax = max(positions(:))+0.1;
    N = size(EigValues,1);
    if j > N
        error('We dont have that many modes!')
    end
    vectors = reshape(EigVectors(:,j),[N/3,3]);
    x = positions(:,1);
    y = positions(:,2);
    z = positions(:,3);
    u = vectors(:,1);
    v = vectors(:,2);
    w = vectors(:,3);
    quiver3(x,y,z,u,v,w,0.3,'r')
    axis equal
    axis([axisMin axisMax axisMin axisMax axisMin axisMax])
end
