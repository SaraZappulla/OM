function [dy] = twobp_perturbed(t,yo,mu,J2,Re)

% This function implements the value of dy = [v dv/dt] = [dr/dt d^2r/dt^2]
% which is the derivated vector of state vector
%
% input 
% t[s]                          time span not necessary for two-body problem
% yo[km km km km/s km/s km/s]   state vector containing initial condition
% mu[km^3/s^2]                  gravitational parameter (depending from the planet/star)
% J2 and Re from astrocontrant 

% output
% dy[km/s km/s km/s km/s^2 km/s^2 km/s^2] derivative of the state vector
% 
% Contributors
% SZ MV
%
% Version
% 1_20/09/2022


ro = yo(1:3); %radius defined from the state vector 
vo = yo(4:6); %velocity defined from the state vector 
vector = [(ro(1)/norm(ro))*((5*(ro(3)^2)/(norm(ro))^2)-1); 
            (ro(2)/norm(ro))*((5*(ro(3)^2)/(norm(ro))^2)-1);
            (ro(3)/norm(ro))*((5*(ro(3)^2)/(norm(ro))^2)-3)];
aj2 = (3/2*J2*mu*Re^2/(norm(ro)^4)).*vector;

dy =[vo; (-mu./(norm(ro).^3).*ro)+aj2]; 

end