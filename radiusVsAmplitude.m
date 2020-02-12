clf
N = 100;
alpha = 1;
[EigValues,EigVectors,positions] = findNormalModes(N,alpha);
radius = sqrt(sum(positions.^2,2));
EigVectors = reshape(EigVectors,N,3,3*N);
amplitude = sqrt(sum(EigVectors.^2,2));
hold on
for i = 1:3*N
    plot(radius,amplitude(:,:,i),'DisplayName',num2str(EigValues(i)));
    
end
legend
hold off
% scatter(radius,amplitude(:,:,1));