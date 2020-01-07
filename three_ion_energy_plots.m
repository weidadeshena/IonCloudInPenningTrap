
%three ions


N=100;

c=1;

%order is x1,y,z1,x2,y2,z2, etc

for i=1:N-1;
    S=i/N;
strength(i)=S;
%here is the energy expression
a=S; %/(1+S);  % axial
b=1-S ; %/(1+S);   %radial


f1=@(y) total_three_1(y,a,b,c);
[y,fval]=fminsearch(f1,[1]);
energy(i,1)= fval;  %  %aligned along z

f2=@(y) total_three_2(y,a,b,c);
[y,energy(i,2)]=fminsearch(f2, [1,1]); 

f2=@(y) total_three_3(y,a,b,c);
[y,energy(i,3)]=fminsearch(f2, [1,1]); 

f2=@(y) total_three_4(y,a,b,c);
[y,energy(i,4)]=fminsearch(f2, [1,1]); 

end

plot( strength, energy(:,[1,2,3,4]))
