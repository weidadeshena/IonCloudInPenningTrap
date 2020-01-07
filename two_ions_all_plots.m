%2 ions


%Here the program plots for all configurations and is told what
%configurations are possible

% S is the strangth patameter - axial strength z
%  S=0 means axial strength = 0
%  S=1 means radial strangth = 0

N=100;
c=1;

%the array x is now as follows
% config 1 x(1) is axial separation
% config 2 x(1) is axial sep, x(2) is centre ion displacement

for i=1:N+1
    S=(i-1)/N;
strength(i)=S;
%here is the energy expression
a=S; %/(1+S);  % axial
b=1-S ; %/(1+S);   %radial

f1=@(x) 2*a*x(1)^2 + (5/2)*c/x(1); %aligned along z
energy1(i,1)= f1(fminsearch(f1, [1]));  %aligned along z

f2=@(x) 2*a*x(1)^2 + 3/2 * b * x(2)^2 + 2*c/((3/2*x(2))^2+x(1)^2)^(1/2)+ c/(2*x(1));
energy2(i,1)= f2(fminsearch(f2, [1,1]));

%f1=@(x) a*x(3)^2 + a*x(6)^2 + c/((x(1)-x(4))^2+(x(2)-x(5))^2+(x(3)-x(6))^2); %aligned along z

%energy(i,2)= f1(fminsearch(f1, [1,-1,.5, -.5, .3, -.3]));  %aligned along z


end

%plot( strength, energy1(:,[1]), strength, energy2(:,[1]))
plot( strength, energy1, strength, energy2)
%plot( strength, energy2(:,[1]))

%plot(ratio, value)

