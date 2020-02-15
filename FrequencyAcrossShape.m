% assume M = 1, Ze^2/(4pi epsilon_0) = 1
clear all
%find radii for N=6,8,30
alpha = [0.1 0.2 0.3 0.4 0.5 0.6 0.85 0.95 1 1.05 1.2 1.67 2 2.5 3.3 5 10];
EigValueData = acrossAlphaFreq(6,alpha);
BreathingModes = sqrt(-min(EigValueData));
plot(alpha,BreathingModes);

function EigValueData = acrossAlphaFreq(N,alpha)
    EigValueData = zeros(3*N,length(alpha));
    for i = 1:length(alpha)
        [EigValues,EigVectors,positions]=findNormalModes(N,alpha(i));
        EigValueData(:,i) = EigValues;
    end
end

function [EigValueData,EigVectorData] = acrossAlphaData(N,alpha)
    EigValueData = zeros(3*N,length(alpha));
    EigVectorData = zeros(3*N,3*N,length(alpha));
    for i = 1:length(alpha)
        [EigValues,EigVectors,wz,w1,r,r_d]=findNormalModes(N,alpha(i));
        EigValueData(:,i) = EigValues;
        EigVectorData(:,:,i) = EigVectors;
    end
end


