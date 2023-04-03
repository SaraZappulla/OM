function [a] = repeating_groundtrack_semimajor_axis(mu,omega_e,m,k)
    omega_e = omega_e*pi/180;
    a = (m^2/((omega_e^2)*(k^2))*mu)^(1/3);
end