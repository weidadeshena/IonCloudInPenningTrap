function energy = total_two_2(y,a,b,c)
% two ions radial
n=2;

energy=0;
x=[y(1),0,0,-y(1),0,0];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);