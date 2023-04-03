%% Earth to Moon direct transfer

clc 
clear all

dep_time1 = date2mjd2000([2009,06,18,05,32,0]); %18 june 2009
arr_time1 = date2mjd2000([2009,06,22,12,0,0]); %22 june 2009
mu1 = 398600; %Earth
mu2 = astroConstants(20); %Moon

%% definition of initial positions

[xP, vP] = ephMoon(arr_time1);
module_xP = norm(xP); %radius module of moon referred to earth

rLEO = 6378; %km of general LEO orbits https://www.space.com/low-earth-orbit
rp = rLEO;  
ra = module_xP;
at = (rp + ra) /2;
et = (ra-rp)/(ra+rp);
pt = at*(1-et^2);
Vat = sqrt(mu1/pt)*(1-et);
Vpt = sqrt(mu1/pt)*(1+et);


module_vP = norm(vP);
vLEO = 0;

deltavp = (Vpt-vLEO)
%deltava = module_vP - vat
tt = (sqrt(at^3/mu1)*pi)/3600/24

%% Insertion in lunar orbit
% respect to moon reference plane the commissioning orbit has these average
% parameters
rmoon = 1737.4; % average moon radius insert website here
ra_com = 216 + rmoon; %km
rp_com = 30 + rmoon; %km
omega_com = 1.5*pi; %rad
e_com = 0.043;
% e = (ra_com-rp_com)/(ra_com+rp_com) verifying correctness data
theta = 0/180*pi; % assumed
a_com = (rp_com + ra_com) /2;
p_com = a_com*(1-e_com^2);

vat = module_vP - Vat;
vr = sqrt(mu2/p_com)*(e_com*sin(theta));
vtheta = sqrt(mu2/p_com)*(1+e_com*cos(theta));
va_com = sqrt(mu2/p_com)*(1-e_com);
v = norm([vr,vtheta]);

deltava = v-vat

vneeded = (567 + 185+133+41)/1000 % from LTI to commissioning

deltava-vneeded

%% orbit representation

%transfer orbit
ToF = tt*24*3600;
tspanT = linspace(0,ToF,10000);
rE = [-6378, 0, 0];
vE = [0 -Vpt 0];
yT = [rE vE];
options= odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[tT,YT] = ode113(@(t,y) twobodyproblem_ode(t,y,mu1), tspanT, yT, options);

figure(1)
hold on
grid on
axis equal
plot3(YT(:,1),YT(:,2),YT(:,3));
plot3(rE(1),rE(2),rE(3),'ro');
plot3(module_xP,0,0,'bo');
earth_sphere

%% From nominal orbit to 
