function energy = total_two_2(y,a,b,c)
% three ions kink
n=3;

energy=0;
x=[y(1),0,y(2),-2*y(1),0,0,y(1),0,-y(2)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);