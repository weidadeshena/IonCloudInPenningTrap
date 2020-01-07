function energy_plots_three_ions

%three ions


N=500;

c=1;

%order is x1,y,z1,x2,y2,z2, etc
strength = zeros(N-1);
energy = zeros(N-1,8);
for i=1:N-1;
    S=i/N;
strength(i)=S;
%here is the energy expression
a=S; %/(1+S);  % axial
b=(1-S) ; %/(1+S);   %radial


f1=@(y) total_four_1(y,a,b,c);
[y,fval]=fminsearch(f1,[1,2]);
energy(i,1)= fval;  %  %aligned along z

f2=@(y) total_four_2(y,a,b,c);
[y,energy(i,2)]=fminsearch(f2, [1,2,1,2]); 

f2=@(y) total_four_3(y,a,b,c);
[y,energy(i,3)]=fminsearch(f2, [1,2,2,1,1]); 

f2=@(y) total_four_4(y,a,b,c);
[y,energy(i,4)]=fminsearch(f2, [1,1]); 

f1=@(y) total_four_5(y,a,b,c);
[y,fval]=fminsearch(f1,[1,2,3]);
energy(i,5)= fval;  %  %aligned along z

f2=@(y) total_four_6(y,a,b,c);
[y,energy(i,6)]=fminsearch(f2, [1,2]); 

f2=@(y) total_four_7(y,a,b,c);
[y,energy(i,7)]=fminsearch(f2, [1]); 

f2=@(y) total_four_8(y,a,b,c);
[y,energy(i,8)]=fminsearch(f2, [1,1,1]); 
end

plot( strength(1:N-1), energy(1:N-1,1:8))
end

function energy = total_four_1(y,a,b,c)
% four ions axial
n=4;

energy=0;
x=[0,0,y(1),0,0,y(2),0,0,-y(2),0,0,-y(1)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);
end

function energy = total_four_2(y,a,b,c)
% four ions kink
n=4;

energy=0;
x=[y(3),0,y(1),y(4),0,y(2),-y(4),0,-y(2),-y(3),0,-y(1)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);
end

function energy = total_four_3(y,a,b,c)
% four ions spiral
n=4;

energy=0;
x=[y(3),y(5),y(1),y(4),0,y(2),-y(4),0,-y(2),-y(3),-y(5),-y(1)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);
end

function energy = total_four_4(y,a,b,c)
% four ions diamond
n=4;

energy=0;
x=[0,0,y(1),y(2),0,0,-y(2),0,0,0,0,-y(1)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);
end

function energy = total_four_5(y,a,b,c)
% four ions offset diamond
n=4;

energy=0;
x=[0,y(3),y(1),y(2),-y(3),0,-y(2),-y(3),0,0,y(3),-y(1)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);
end

function energy = total_four_6(y,a,b,c)
% four ions offset square
n=4;

energy=0;
x=[y(1),0,y(2),-y(1),0,y(2),0,y(1),-y(2),0,-y(1),-y(2)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);
end

function energy = total_four_7(y,a,b,c)
% four ions square
n=4;

energy=0;
x=[y(1),0,0,-y(1),0,0,0,y(1),0,0,-y(1),0];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);
end

function energy = total_four_8(y,a,b,c)
% four ions tetrahedron
n=4;

energy=0;
x=[0,0,3*y(1),0,2*y(3),-y(1),y(2),-y(3),-y(1),-y(2),-y(3),-y(1)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);
end