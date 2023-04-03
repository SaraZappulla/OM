function [dy] = twobp(t, yo,mu)

% function that implement the value of dy = [v dv/dt] = [dr/dt d^2r/dt^2]
%
% input 
% t[s]     time span not necessary for two-body problem
% yo[km km km km/s km/s km/s] state vector containing initial condition
% mu[km^3/s^2] gravitational parameter (depending from the planet/star)
% 
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

dy =[vo; -mu./(norm(ro).^3).*ro]; 

end
