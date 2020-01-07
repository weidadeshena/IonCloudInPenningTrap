alpha = 3;
N = 4;
G = Greens(alpha);
omega_z_squared = G;
omega_1_squared = (1-G)/2;
final_positions = crystal_graphs_energy(N,omega_z_squared,omega_1_squared);

 mm,
function G = Greens(alpha)
    if alpha < 1
        beta = sqrt(1-alpha^2);
        G = 1/beta^2 - alpha/beta^3*atan(beta/alpha);
    elseif alpha > 1
        beta = sqrt(alpha^2-1);
        G = -1/beta^2 + alpha/(2*beta^3)*log(abs((alpha+beta)/(alpha-beta)));
    end
end
    