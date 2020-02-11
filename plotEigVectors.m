clear all
clf
N = 60;
alpha = 1;
repeat = 50;

data = zeros(3*N,repeat);
for i = 1:repeat
    [EigValues,EigVectors,positions] = findNormalModes(N,alpha);
    data(:,i) = EigValues;
end


% clf
% hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
%              'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
%              'XLim',[-1.01 0.01],...               %#   set the x axis limit,
%              'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
%              'Color','none');               %#   and don't use a background color
% plot(EigValues,0,'bx');

histogram(data,50)


% folder = [N,'IonModes_alpha',num2str(alpha)];
% if 7 ~= exist(folder,'dir')
%     mkdir(folder);
% end
% for i = 1:3*N
%     fig = plotMode(i);
% %     temp = ['Mode',num2str(i),'freq',num2str(EigValues(i)),'.fig'];
% %     temp = fullfile(folder,temp);
% %     saveas(fig,temp)
% end

% plotMode(963)

function fig = plotMode(j)
    global EigVectors positions N;
% global EigValues EigVectors positions;
    axisMin = min(positions(:))-0.3;
    axisMax = max(positions(:))+0.3;
    if j > 3*N
        error('We dont have that many modes!')
    end
    vectors = reshape(EigVectors(:,j),[N,3]);
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
