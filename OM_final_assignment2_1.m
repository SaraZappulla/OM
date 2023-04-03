%% Project model:

% 1. nominal orbit representation
% 2. nominal orbit groud track for 1 orbit, 1 day and 10 days, unperturbed
% 3. modify a to repeat ground track for unperturbed using ratio given
% 4. plot groundtrack nominal and modified adding perturbations
% 5. propagation of orbit with perturbation: expressed in cartesian
% coordinate and keplerian elements
% 6. plot all the elements both propagation methods, in terms of error,
% computation time...
% 7. represent orbit evolution
% 8. filtering of high frequenies and plotting of both
% 9. select an object and download the orbital elements - propagate its
% orbit using our model - compare real and our model

% Points 1-4 are in OM_final_assignmenta_2_1
% Points 5-9 are in OM_final_assignmenta_2_2 project assignment data

%% orbit data
kep = [22594, 0.6847, 63.6246, 0,0,0];

%SRP data
p_sr = 4.5e-6; %at 1 AU, N/m^2
cr = 1.0; %reflection coefficient
Am_ratio = 5; %surface to mass ratio relative to the sun, m^2/kg 

%writing the functions
obl = astroConstants(8); %Earth's obliquity
mu = astroConstants(13); %earth
[r_spacecraft_earth,v_spacecraft_earth] = kepl_to_car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),mu);
r_spacecraft_sun = r_spacecraft_earth + [cos(obl);sin(obl);0]*astroConstants(2);
acc_SRP = -p_sr*cr*Am_ratio.*r_spacecraft_sun./norm(r_spacecraft_sun);

%degrees to radians conversion, to be done beforehand
kep_rad = [kep(1) kep(2) 0 0 0 0];
kep_rad(3:6) = kep(3:6)*pi/180;

Tperiod = 2*pi * sqrt(kep(1)^3/mu);
tspan = linspace(0,Tperiod*100,1e5);

acc_perturbed = a_per_complete(tspan,kep_rad,acc_SRP); %RSW ref. frame
%ds_per = EOM(tspan,kep_rad,mu,acc_perturbed);

% s = state in this case written keplerian parameters
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,kep_propagate] = ode113(@(t,s) (EOM(t,s,mu,a_per(t,s))'), tspan, kep_rad, options);

%rotation of acc_perturbed from RSW to carthesian ref. frame

%% 1. Nominal orbit representation


yo = [r_spacecraft_earth,v_spacecraft_earth];
tspan = linspace(0,Tperiod,10000); %%time span where the state vector is defined
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y] = ode113(@(t,y) twobodyproblem_ode(t,y,mu), tspan, yo, options);

figure(1)
earth_sphere;
axis equal
grid on
hold on;
plot3(Y(:,1),Y(:,2),Y(:,3),'k-');
plot3(yo(1),yo(2),yo(3),'ro');
legend('Earth','Nominal orbit','Spacecraft initial position')
title('Nominal orbit [km]')

%% 2. Nominal orbit groud track for 10 days, unperturbed

%calculation of the orbital period
time_vector = (linspace(0,24*3600*10,100000))'; 

%groundtrack function
[alpha,delta,lon,lat] = groundtrack(mu,thetag_zero,time_vector,omega_e,kep(1),kep(2),kep(3),kep(4),kep(5),kep(6));
lon180 = wrapTo180(lon);
lat180 = wrapTo180(lat);

%singularity plotting conditions on latitude and longitude
lat(end+1) = inf;
lon(end+1) = inf;
figure(1)
xlin = [-180 180];
ylin = [90 -90]; %do 90 -90 because otherwise we get axis upside down
I = imread("Earth for groundtracks.jpg","jpg");
h = image('Xdata',xlin,'Ydata',ylin,'Cdata',I);
hold on
plot(lon180,lat180,'LineStyle','none','Marker','.','Color','p','MarkerSize',2)
% title('Ground track 10 day, unperturbed problem')
% xlabel('Longitude [°]')
% ylabel('Latitude [°]')

%PA PUT THE GRAPH SAME DIMENSION NOT 100 OR 200

%% 2. Nominal orbit groud track for 1 day, unperturbed
%calculation of the orbital period
time_vector = (linspace(0,24*3600,100000))'; 

%groundtrack function
[alpha,delta,lon,lat] = groundtrack(mu,thetag_zero,time_vector,omega_e,kep(1),kep(2),kep(3),kep(4),kep(5),kep(6));
lon180 = wrapTo180(lon);
lat180 = wrapTo180(lat);

%singularity plotting conditions on latitude and longitude
lat(end+1) = inf;
lon(end+1) = inf;
figure(1)
xlin = [-180 180];
ylin = [90 -90]; %do 90 -90 because otherwise we get axis upside down
% I = imread("Earth for groundtracks.jpg","jpg");
% h = image('Xdata',xlin,'Ydata',ylin,'Cdata',I);
hold on
plot(lon180,lat180,'LineStyle','none','Marker','.','Color','r','MarkerSize',2)
% title('Ground track 1 day, unperturbed problem')
% xlabel('Longitude [°]')
% ylabel('Latitude [°]')

%PA PUT THE GRAPH SAME DIMENSION NOT 100 OR 200
%% 2. Nominal orbit groud track for 1 orbit, unperturbed

thetag_zero = 0; %degrees
time_vector = (linspace(0,80000,9000))'; 
omega_e = 15.04/3600; %deg/s, depends on time_vector spacing

%calculation of the orbital period
time_vector = (linspace(0,Tperiod,100000))'; 

%groundtrack function
[alpha,delta,lon,lat] = groundtrack(mu,thetag_zero,time_vector,omega_e,kep(1),kep(2),kep(3),kep(4),kep(5),kep(6));
lon180 = wrapTo180(lon);
lat180 = wrapTo180(lat);

%singularity plotting conditions on latitude and longitude
lat(end+1) = inf;
lon(end+1) = inf;
figure(1)
% xlin = [-180 180];
% ylin = [90 -90]; %do 90 -90 because otherwise we get axis upside down
% I = imread("Earth for groundtracks.jpg","jpg");
% h = image('Xdata',xlin,'Ydata',ylin,'Cdata',I);
hold on
plot(lon180,lat180,'LineStyle','none','Marker','.','Color','g','MarkerSize',2)
% title('Ground track 1 orbit, unperturbed problem')
% xlabel('Longitude [°]')
% ylabel('Latitude [°]')

% UNIFY AND ADD LEGEND

%% 3. Modify a to repeat ground track for unperturbed using ratio given

thetag_zero = 0; %degrees
time_vector = (linspace(0,80000,9000))'; 
omega_e = 15.04/3600; %deg/s, depends on time_vector spacing

%te = 365* %earth period
%Tearth = insert here the spacecraft ratio that comes from 13:5 respect to
%earth
%calculation of the orbital period
% time_vector = (linspace(0,Tearth,100000))'; 
% 
% %groundtrack function
% [alpha,delta,lon,lat] = groundtrack(mu,thetag_zero,time_vector,omega_e,kep(1),kep(2),kep(3),kep(4),kep(5),kep(6));
% lon180 = wrapTo180(lon);
% lat180 = wrapTo180(lat);
% 
% %singularity plotting conditions on latitude and longitude
% lat(end+1) = inf;
% lon(end+1) = inf;
% figure(1)
% xlin = [-180 180];
% ylin = [90 -90]; %do 90 -90 because otherwise we get axis upside down
% I = imread("Earth for groundtracks.jpg","jpg");
% h = image('Xdata',xlin,'Ydata',ylin,'Cdata',I);
% hold on
% plot(lon180,lat180,'LineStyle','none','Marker','.','Color','g','MarkerSize',2)
% title('Ground track 1 orbit, unperturbed problem')
% xlabel('Longitude [°]')
% ylabel('Latitude [°]')

%% 4. plot groundtrack nominal and modified adding perturbations

%degrees to radians conversion, to be done beforehand
kep_rad = [kep(1) kep(2) 0 0 0 0];
kep_rad(3:6) = kep(3:6)*pi/180;

Tperiod = 2*pi * sqrt(kep(1)^3/mu);
tspan = linspace(0,Tperiod*100,1e5);

acc_perturbed = a_per_complete(tspan,kep_rad,acc_SRP); %RSW ref. frame
%ds_per = EOM(tspan,kep_rad,mu,acc_perturbed);

% s = state in this case written keplerian parameters
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,kep_propagate] = ode113(@(t,s) (EOM(t,s,mu,a_per_complete(t,s,acc_SRP))'), tspan, kep_rad, options);
%calculation of the orbital period
time_vector = (linspace(0,10000,100000))'; 

thetag_zero = 0; %degrees
omega_e = 1000/3600; %deg/s, depends on time_vector spacing

%groundtrack function
%kep_propagate(3:6,1) = kep_propagate(3:6,1)*180/pi;
[alpha,delta,lon,lat] = groundtrack(mu,thetag_zero,time_vector,omega_e,kep(1),kep(2),kep(3),kep(4),kep(5),kep(6));
lon180 = wrapTo180(lon);
lat180 = wrapTo180(lat);

[alpha,delta,lon_per,lat_per] = groundtrack(mu,thetag_zero,time_vector,omega_e,...
    kep_propagate(end,1),kep_propagate(end,2),kep_propagate(end,3),kep_propagate(end,4),kep_propagate(end,5),kep_propagate(end,6));
lon180_per = wrapTo180(lon_per);
lat180_per = wrapTo180(lat_per);



%singularity plotting conditions on latitude and longitude
lat_per(end+1) = inf;
lon_per(end+1) = inf;
lat(end+1) = inf;
lon(end+1) = inf;
figure(4)
xlin = [-180 180];
ylin = [90 -90]; %do 90 -90 because otherwise we get axis upside down
I = imread("Earth for groundtracks.jpg","jpg");
h = image('Xdata',xlin,'Ydata',ylin,'Cdata',I);
hold on
plot(lon180,lat180,'LineStyle','none','Marker','.','Color','g','MarkerSize',2)
plot(lon180_per,lat180_per,'LineStyle','none','Marker','.','Color','r','MarkerSize',2)
title('Ground track')
xlabel('Longitude [°]')
ylabel('Latitude [°]')
legend('unperturbed','perturbed')
% GRAPH TO FIX
