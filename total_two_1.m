function energy = total_two_1(y,a,b,c)
% two ions
n=2;
energy=0;
x=[0,0,y(1),0,0,-y(1)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);