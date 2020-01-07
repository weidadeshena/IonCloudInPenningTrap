function energy = total_two_2(y,a,b,c)
% three ions axial
n=3;

energy=0;
x=[0,0,y(1),0,0,0,0,0,-y(1)];
energy = coulomb_energy(n,x,c) + trap_energy(n,x,a,b);