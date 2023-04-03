%% Exercise 2

%Data definition
mu_planet = astroConstants(13); %Earth case
mu_sun = astroConstants(4);
AU = astroConstants(2);
V_neg = [31.5;5.2;0]; %km/s
V_pos = [36;0;0]; %km/s
r_earth = [0;-1;0]*AU; %Earth's position, km
R_earth = 6371; %Earth's radius, km
h_atm = 100; %Earth's atmosphere altitude, km 
%computing Vplanet
v_planet_modulus = sqrt(mu_sun/norm(r_earth));
v_planet = v_planet_modulus.*[1;0;0]; %oriented in the X direction

%calculating approach and departure velocities
v_inf_neg = V_neg-v_planet;
v_inf_pos = V_pos-v_planet;

%calculating the turning angle
v_neg = norm(v_inf_neg);
v_pos = norm(v_inf_pos);

turning_angle = acos(dot(v_inf_neg,v_inf_pos)/(v_pos*v_neg));

%calculating r_peri
turning_angle_function = @(r_peri) (2*asin(1/(1+r_peri*(v_neg^2)/mu_planet)))/2+...
    (2*asin(1/(1+r_peri*(v_pos^2)/mu_planet)))/2-turning_angle;

%definition of an initial guess
r_peri_guess = 7000; %Earth's radius + random height, km
r_peri_fzero = fzero(turning_angle_function,r_peri_guess);
if r_peri_fzero < (R_earth+h_atm)
    fprintf('error, too low pericenter radius');
end

%pericenter velocities
e_in = 1+r_peri_fzero*(v_neg^2)/mu_planet; %incoming arc eccentricity
e_out = 1+r_peri_fzero*(v_pos^2)/mu_planet; %outgoing arc eccentricity
a_in = -mu_planet/(v_neg^2); %incoming arc semimajor axis, km
a_out = -mu_planet/(v_pos^2); %outgoing arc semimajor axis, km
p_in = a_in*(1-e_in^2); %incoming arc semilatus rectum, km
p_out = a_out*(1-e_out^2); 
v_peri_in = sqrt(mu_planet/p_in)*(1+e_in);
v_peri_out = sqrt(mu_planet/p_out)*(1+e_out);
deltaV_peri = abs(v_peri_out-v_peri_in);

%calculating deltaV of the entire flyby
deltaV_flyby = norm(V_pos-V_neg);

%arcs plotting, planetocentric
r_peri = r_peri_fzero.*[1;0;0]; 
yo = [r_peri v_peri_in.*[0;1;0]]; %arrival state vector (pericenter)
y1 = [r_peri v_peri_out.*[0;1;0]]; %departing state vector (pericenter)
%T = 10000000; 
tspano = linspace(10000,0,10000); %time span where the state vector is defined
tspan1 = linspace(0,10000,10000);
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_planet), tspan1, y1, options);
[T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_planet), tspano, yo, options);
Yo = Yo./R_earth;
Y1 = Y1./R_earth;
figure(5);
grid on
hold on
axis equal
plot3(Yo(:,1),Yo(:,2),Yo(:,3),'b-');
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');

%earth plot with Earth's radius scale-adjusted dimensions
earth_sphere_earth_radius_scaled

%asymptotes plotting - issue on the angular coefficients
turning_angle_incoming = asin(1/e_in);
turning_angle_outcoming = asin(1/e_out);

x_incoming = -linspace(0,15*R_earth,10000)*cos(pi/2-turning_angle_incoming)/R_earth-a_in/R_earth+1;
y_incoming = -linspace(0,15*R_earth,10000)*sin(pi/2-turning_angle_incoming)/R_earth;
z_incoming = -linspace(0,15*R_earth,10000)*v_inf_neg(3)/R_earth;
plot3(x_incoming,y_incoming,z_incoming)

x_outgoing = linspace(0,15*R_earth,10000)*cos(pi/2+turning_angle_outcoming)/R_earth-a_out/R_earth+1;
y_outgoing = linspace(0,15*R_earth,10000)*sin(pi/2+turning_angle_outcoming)/R_earth;
z_outgoing = linspace(0,15*R_earth,10000).*v_inf_pos(3)./R_earth;
plot3(x_outgoing,y_outgoing,z_outgoing)
