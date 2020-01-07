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
% config 3 x(1) is radial, x(2) is axial
%config 4 triangle x(1) is radius



for i=1:N-1
    S=(i)/N;
strength(i)=S;


if S>=1
    S=.999
end

a=S;   % axial
b=1-S ;  %radial


    

f1=@(x) 2*a*x(1)^2 + (5/2)*c/x(1); %aligned along z
energy1(i,1)= f1(fminsearch(f1, [1]));  %aligned along z


f2=@(x) 2*a*x(1)^2 + 3/2 * b * x(2)^2 + 2*c/((3/2*x(2))^2+x(1)^2)^(1/2)+ c/(2*x(1));
energy2(i,1)= f2(fminsearch(f2, [1,1]));

f3=@(x) 2*b*x(1)^2 + 3/2 * a * x(2)^2 + 2*c/((3/2*x(2))^2+x(1)^2)^(1/2)+ c/(2*x(1));
energy3(i,1)= f3(fminsearch(f3, [1,1]));


f4=@(x) 3*b*x(1) +3*c/(3^(1/2)*x(1));
energy4(i,1)= f4(fminsearch(f4, [1]));



end

plot( strength, energy1, strength, energy2, strength, energy3, strength, energy3)
ylim([0 5])
xlim([0,1])




