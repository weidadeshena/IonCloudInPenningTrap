% clear all
clf
global EigValues EigVectors positions N;
% load NormalModes EigValues EigVectors positions
load NormalModes_nonspherical EigValues EigVectors positions


N = size(EigValues,1);

% folder = [num2str(N/3),'IonModes'];
folder = [num2str(N/3),'IonModes_alpha0.1'];
if 7 ~= exist(folder,'dir')
    mkdir(folder);
end
for i = 1:N
    fig = plotMode(i);
    temp = ['Mode',num2str(i),'freq',num2str(EigValues(i)),'.fig'];
    temp = fullfile(folder,temp);
    saveas(fig,temp)
end

% plotMode(963)

function fig = plotMode(j)
    global EigVectors positions N;
% global EigValues EigVectors positions;
    axisMin = min(positions(:))-0.3;
    axisMax = max(positions(:))+0.3;
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
    fig = plot3(x,y,z,'b*','MarkerSize',4);
    hold on
    quiver3(x,y,z,u,v,w,0.3,'r');
    hold off
    axis equal
    axis([axisMin axisMax axisMin axisMax axisMin axisMax])
end
