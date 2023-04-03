function[alpha,delta,lon,lat] = groundtrack(mu,thetag_zero,time_vector,omega_e,a,e,i,OM,om,theta,rr,vv)

if(nargin==10) 
    %calculation of initial carthesian conditions
    [rr,vv] = kepl_to_car(a,e,i,OM,om,theta,mu);
end
if(nargin==6)
    rr=a;
    vv=e;
    %calculation of keplerian elements (unnecessary)
    [a,e,i,OM,om,theta] = car_to_kepl(rr,vv,mu);
end

%Initial integration conditions
s0 = [rr;vv];

%Numerical integration error, upper limit
options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[~, state] = ode45(@(t,s) twobodyproblem_ode(t,s,mu),time_vector,s0,options);

rt = state(:,1:3); 
delta = asin(rt(:,3)./[(vecnorm(rt'))]'); %time-changing vector, in radians
alpha = atan2(rt(:,2),rt(:,1)); %radians

%Radians to degrees conversion, necessary because omega_e and thetag_zero
%are provided in degrees
alpha = alpha.*(180/pi); 
delta = delta.*(180/pi);

theta_g = thetag_zero + omega_e.*time_vector;

lon = alpha - theta_g;
lat = delta;

end
