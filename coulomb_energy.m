function energy=coulomb_energy(n,x,c)
energy=0;
for i=1:n-1
    for j=i+1:n
        xx=x(3*i-2)-x(3*j-2);
        yy=x(3*i-1)-x(3*j-1);
        zz=x(3*i)-x(3*j);
        d=sqrt(xx^2+yy^2+zz^2);
        energy=energy+c/d;
    end
end

        