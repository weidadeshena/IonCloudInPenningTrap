% the amplitude of a specific normal modes on a ion plotted against the 
% radial/axial position of the ion
N = 80;
alpha = 1;
[EigValues,EigVectors,positions] = findNormalModes(N,alpha);
EigVectors = reshape(EigVectors,N,3,3*N);
lowerRange = -1;
upperRange = -0.98;

folder = [num2str(N),'Ion_alpha',num2str(alpha),'range',num2str(lowerRange),num2str(upperRange)];
plotRadialAxial(positions,EigValues,EigVectors,lowerRange,upperRange,folder)


function plotAmplitude(positions,EigVectors)
    radius = sqrt(sum(positions.^2,2));
    amplitude = sqrt(sum(EigVectors.^2,2));
    hold on

    for i = 1:3*N
        plot(radius,amplitude(:,:,i),'-x','DisplayName',num2str(EigValues(i)));

    end
    legend
    hold off
    % scatter(radius,amplitude(:,:,1));
end



% check if a mode is radially compressed and axially expand
function plotRadialAxial(positions,EigValues,EigVectors,lowerRange,UpperRange,folder)
    [sortedEigValues,orderIndex] = sort(EigValues);
    sortedEigVectors = EigVectors(:,:,orderIndex);
    lowerIndex = find(sortedEigValues>lowerRange,1);
    upperIndex = length(EigValues) - find(flip(sortedEigValues)<UpperRange,1) + 1;
    EigValuesRange = sortedEigValues(lowerIndex:upperIndex);
    radial = sqrt(sum(positions(:,[1 2]).^2,2));
    radialAmplitude = sqrt(sum(sortedEigVectors(:,[1,2],lowerIndex:upperIndex).^2,2));
    axial = positions(:,3);
    axialAmplitude = sortedEigVectors(:,3,lowerIndex:upperIndex);
    maxAmplitude = max([max(radialAmplitude) max(axialAmplitude)]);
    if 7 ~= exist(folder,'dir')
        mkdir(folder);
    end
    
    for i = 1:(upperIndex-lowerIndex+1)
        fig = figure;
        subplot(1,2,1)
        scatter(radial,radialAmplitude(:,:,i),'DisplayName',num2str(EigValuesRange(i)));
        legend
        ylim([-maxAmplitude(i)-0.05 maxAmplitude(i)+0.05]);
        subplot(1,2,2)
        scatter(axial,axialAmplitude(:,:,i),'DisplayName',num2str(EigValuesRange(i)));
        legend
        ylim([-maxAmplitude(i)-0.05 maxAmplitude(i)+0.05]);
        temp = ['freq',num2str(EigValuesRange(i)),'.fig'];
        temp = fullfile(folder,temp);
        saveas(fig,temp)
    end
    
    
end

