function [omega_z_squared,omega_1_squared] = Greens(alpha)
% calculate the Greens function from Y.Kiwamoto 2004
    if alpha < 1
        beta = sqrt(1-alpha^2);
        G = 1/beta^2 - alpha/beta^3*atan(beta/alpha);
    elseif alpha > 1
        beta = sqrt(alpha^2-1);
        G = -1/beta^2 + alpha/(2*beta^3)*log(abs((alpha+beta)/(alpha-beta)));
    else
        G = 1/3;
    end
    omega_z_squared = G;
    omega_1_squared = (1-G)/2;
end