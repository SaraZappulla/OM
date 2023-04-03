function ds = twobodyproblem_ode(~,s,mu)
%ODE (Cauchy problem) solver for the restricted two body problem
%
%INPUTS
%t: time {1x1} (~ for autonomous problems)
%s: state vector [r;v] {6x1}
%mu: gravitational parameter of Earth {1x1}
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

%output definition
ds=[v;-mu/r_mod^3*r(1);-mu/r_mod^3*r(2);-mu/r_mod^3*r(3)];
end
