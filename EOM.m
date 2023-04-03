function [kep_dot] = EOM(t,kep,mu,acc_per)

%acc_per = ar, as, aw
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);
ar = acc_per(1,1);
as = acc_per(2,1);
aw = acc_per(3,1);

% evertything needed
p = a*(1-e^2);
r = p/(1+e*cos(theta));
v = sqrt(2*mu/r-mu/a);
h = sqrt(p*mu);
b = a*sqrt(1-e^2);
n = sqrt(mu/(a^3));

% effective equations
a_dot = 2*a^2/h*(e*sin(theta)*ar + p/r*as);
e_dot = 1/h * (p*sin(theta)*ar + ((p+r)*cos(theta)+r*e)*as);
i_dot = r*cos(theta+om)/h*aw;
OM_dot = r*sin(theta+om)/(h*sin(i))*aw;
om_dot = 1/(h*e) * (-p*cos(theta)*ar + (p+r)*sin(theta)*as)-r*sin(theta+om)*cos(i)/(h*sin(i))*aw;
theta_dot = h/(r^2) + 1/(h*e) * (p*cos(theta)*ar - (p+r)*sin(theta)*as);

kep_dot = [a_dot e_dot i_dot OM_dot om_dot theta_dot];
%[rr_per vv_per]=kepl_to_car(a_dot,e_dot,i_dot*180/pi,OM_dot*180/pi,om_dot*180/pi,theta_dot*180/pi,mu);
%ds_per = [rr_per;vv_per];
end
