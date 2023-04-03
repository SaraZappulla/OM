function ds = twobodyproblem_perturbed_ode(~,s,mu,re,jj,a_SRP)
%ODE (Cauchy problem) solver for the restricted two body problem with the
%J2 zonal perturbation element
%
%INPUTS
%t: time {1x1} (~ for autonomous problems)
%s: state vector [r;v] {6x1}
%mu: gravitational parameter of Earth {1x1}
%re: equatorial radius of Earth {1x1}
%jj: constant coefficient {1x1}
%
%OUTPUT
%ds: derivative of the state vector in respect to time [rdot,vdot] {6x1}
%
%-----------------------------------------------------------------------

%definition of s as a column vector
r=s(1:3);
v=s(4:6);

%definition of the lenght of vector r (2-norm)
r_mod=norm(r);

%output definition, carthesian coordinates
ds=[v;a_SRP(1)-mu/r_mod^3*r(1)+((3/2*jj*mu*(re^2))/(r_mod^4))*(r(1)/r_mod)*(5*((r(3)/(r_mod))^2)-1);...
    a_SRP(2)-mu/r_mod^3*r(2)+((3/2*jj*mu*(re^2))/(r_mod^4))*(r(2)/r_mod)*(5*((r(3)/(r_mod))^2)-1);...
    a_SRP(3)-mu/r_mod^3*r(3)+((3/2*jj*mu*(re^2))/(r_mod^4))*(r(3)/r_mod)*(5*((r(3)/(r_mod))^2)-3)];

end

