%4 ions


%Here the program chooses which configuration for three ions without being
%told which is the lowest energy

% S is the strangth patameter - axial strength z
%  S=0 means axial strength = 0
%  S=1 means radial strangth = 0

N=100;
c=1;

%order is x1,y,z1,x2,y2,z2, etc

for i=1:N+1
    S=(i-1)/N;
strength(i)=S;
%here is the energy expression
a=S; %/(1+S);  % axial
b=1-S ; %/(1+S);   %radial


f1=@(x) a*(x(3)^2 + x(6)^2 + x(9)^2) + x(12)^2 + b*(1.01*x(1)^2+x(2)^2+x(4)^2+x(5)^2+x(7)^2+x(8)^2+ x(10)^2+x(11)^2)+  c/sqrt(((x(1)-x(4))^2+(x(2)-x(5))^2+(x(3)-x(6))^2))+  c/sqrt(((x(1)-x(7))^2+(x(2)-x(8))^2+(x(3)-x(9))^2)) + c/sqrt(((x(7)-x(4))^2+(x(8)-x(5))^2+(x(9)-x(6))^2)); %aligned along z
energy(i,1)= f1(fminsearch(f1, [1,-1,.5, -.5, .3, -.3, .8, .6, .4, 1.3, 2.4, 1.6]));  %aligned along z

%f1=@(x) a*x(3)^2 + a*x(6)^2 + c/((x(1)-x(4))^2+(x(2)-x(5))^2+(x(3)-x(6))^2); %aligned along z

%energy(i,2)= f1(fminsearch(f1, [1,-1,.5, -.5, .3, -.3]));  %aligned along z


end

plot( strength, energy(:,[1]))
%plot(ratio, value)

