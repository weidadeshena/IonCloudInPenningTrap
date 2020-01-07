function energy=trap_energy(n,x,a,b)

% n = no of ions
% x = positions x1,y1,z1,x2,y2,z2 etc
% a = radial strength
% b = axial strength

% returns total trap potential energy

energy=0;
for i=1:n
    energy=energy+a*x(3*i)^2+b*(x(3*i-1)^2+x(3*i-2)^2);
end
