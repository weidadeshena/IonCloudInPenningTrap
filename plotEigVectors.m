% create a folder with name convention as "[number of
% ion]IonModes_alpha[value of alpha]", find the eigenvectors using
% findNormalModes.m and plot it out and store it in the folder in .fig
% format

clf
N = 50;
alpha = 0.1;

% repeat = 50;
% data = zeros(3*N,repeat);
% for i = 1:repeat
%     [EigValues,EigVectors,positions] = findNormalModes(N,alpha);
%     data(:,i) = EigValues;
% end
% N_list = [2,3,4,5,6,7,8,10,15,20,30,40,50,75,100];
% N_list = [2,3]
% N_size = size(N_list,2);
% tiledlayout(N_size,1)
% for i = 1:N_size
%     [EigValues,EigVectors,positions] = findNormalModes(N_list(i),alpha);
%     data = EigValues(EigValues < 0);
%     nexttile
%     hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
%                  'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
%                  'XLim',[-0.01 1.01],...               %#   set the x axis limit,
%                  'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
%                  'Color','none');               %#   and don't use a background color
%     plot(data,0,'bx');
%     xlabel('\omega^2')
%     temp = ['N = ',N_list(i)];
%     ylabel('temp')
% end
% histogram(data,50)


folder = [num2str(N),'IonModes_alpha',num2str(alpha)];
if 7 ~= exist(folder,'dir')
    mkdir(folder);
end
[EigValues,EigVectors,positions] = findNormalModes(N,alpha);
for i = 1:3*N
    fig = plotMode(i,EigVectors,positions,N);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    set(gca,'Box','on')
    temp = ['Mode',num2str(i),'freq',num2str(EigValues(i)),'.fig'];
    temp = fullfile(folder,temp);
    saveas(fig,temp)
end

% plotMode(963)

function fig = plotMode(j,EigVectors,positions,N)
    axisMin = min(positions(:))-0.5;
    axisMax = max(positions(:))+0.5;
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
    quiver3(x,y,z,u,v,w,0.5,'r');
    
    hold off
    axis equal
    axis([axisMin axisMax axisMin axisMax axisMin axisMax])
    grid on
end
