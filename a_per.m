function [acc_perturbed] = a_per(t,kep)

%NOTE: if the planet is not earth change them
j2 = astroConstants(9);
mu = astroConstants(13);
R = astroConstants(23);

i = kep(3);
omega = kep(5);
theta = kep(6);

%perifocal_to_car frame
[r,V] = kepl_to_car(kep(1),kep(2),kep(3)*180/pi,kep(4)*180/pi,kep(5)*180/pi,kep(6)*180/pi,mu);

vector = [1 - 3*sin(i)^2*sin(theta+omega)^2; ...
            sin(i)^2*sin(2*(theta+omega));...
            sin(2*i)*sin(theta+omega)];
%RSW reference frame
r = norm(r);
constant_term = -3/2 * (j2.*mu.*R.^2) ./(r.^4);
acc_perturbed = constant_term .* vector;
end