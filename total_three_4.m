function energy = total_two_3(y,a,b,c)
% three ions radial triangle
n=3;
energy=0;
x=[y(1),y(2),0,-y(1),y(2),0,0,-2*y(2),0];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);