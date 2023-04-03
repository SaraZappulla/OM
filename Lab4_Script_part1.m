%To do at the end: defining T, improving figures

%% Exercise 1a

%Data definition
mu_planet = astroConstants(13); %Earth case
mu_sun = astroConstants(4);
AU = astroConstants(2);
v_inf_neg = [15.1;0;0]; %km/s
impact_par_modulus = 9200; %km
r_earth = [1;0;0]*AU; %Earth's position, km

%hyperbolic trajectory characterization
v_inf = norm(v_inf_neg); %modulus of v_inf_neg, km/s
a = -mu_planet/(v_inf^2); %semimajor axis, km
turning_angle = 2*atan(-a/impact_par_modulus); %Rad
e = 1/sin(turning_angle/2); %eccentricity
r_peri_modulus = a*(1-e); %pericenter radius, km

%computing deltaV modulus
deltaV_modulus = 2*v_inf*sin(turning_angle/2); %km/s

%% Passage in front of the planet (case 1)

%assumption: u = [0;0;1]
u = [0;0;1];

%rotating the approach velocity vector
v_inf_pos = rodrigues(v_inf_neg,u,-turning_angle); %negative angle, leading

v_planet_modulus = sqrt(mu_sun/norm(r_earth));
v_planet = v_planet_modulus.*[0;1;0]; %oriented in the Y direction

%computing capital V's (vectors in the heliocentric frame)
V_pos = v_planet+v_inf_pos;
V_neg = v_planet+v_inf_neg;

%PLOTTING WORKS, HAS TO BE IMPROVED

%We have to use a negative timespan on the planet approaching hyperbola
%(positive time going backwards); the entire timespan is unknown

%conditions for heliocentric arcs plotting
r_peri = r_peri_modulus.*[0;1;0];
yo = [r_earth+r_peri V_neg]; %arrival state vector
y1 = [r_earth+r_peri V_pos]; %departing state vector
T = 10000000; 
tspano = linspace(20000000,0,10000); %time span where the state vector is defined
tspan1 = linspace(0,20000000,10000);
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_sun), tspan1, y1, options);
[T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_sun), tspano, yo, options);
Y1 = Y1./AU;
Yo = Yo./AU;
figure(1)
grid on
hold on;
plot3(Yo(:,1),Yo(:,2),Yo(:,3),'b-');
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(0,0,0,'yo');

%% Passage behind the planet (case 2)

%assumption: u = [0;0;1]
u = [0;0;1];

%rotating the approach velocity vector
v_inf_pos = rodrigues(v_inf_neg,u,turning_angle); %positive angle, trailing

v_planet_modulus = sqrt(mu_sun/norm(r_earth));
v_planet = v_planet_modulus.*[0;1;0]; %oriented in the Y direction

%computing capital V's (vectors in the heliocentric frame)
V_pos = v_planet+v_inf_pos;
V_neg = v_planet+v_inf_neg;

%PLOTTING WORKS, HAS TO BE IMPROVED

%We have to use a negative timespan on the planet approaching hyperbola
%(positive time going backwards); the entire timespan is unknown

%conditions for heliocentric arcs plotting
r_peri = r_peri_modulus.*[0;1;0];
yo = [r_earth+r_peri V_neg]; %arrival state vector
y1 = [r_earth+r_peri V_pos]; %departing state vector

%defining the integration time T [to be done, maybe at the end of the lab]

tspano = linspace(20000000,0,10000); %time span where the state vector is defined
tspan1 = linspace(0,20000000,10000);
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_sun), tspan1, y1, options);
[T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_sun), tspano, yo, options);
Y1 = Y1./AU;
Yo = Yo./AU;
figure(2)
grid on
hold on;
plot3(Yo(:,1),Yo(:,2),Yo(:,3),'b-');
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(0,0,0,'yo','LineWidth',2);
%% Passage under the planet (case 3)

%assumption: u = [0;1;0]
u = [0;1;0];

%rotating the approach velocity vector
v_inf_pos = rodrigues(v_inf_neg,u,-turning_angle); %negative angle

v_planet_modulus = sqrt(mu_sun/norm(r_earth));
v_planet = v_planet_modulus.*[0;1;0]; %oriented in the Y direction

%computing capital V's (vectors in the heliocentric frame)
V_pos = v_planet+v_inf_pos;
V_neg = v_planet+v_inf_neg;

%PLOTTING WORKS, HAS TO BE IMPROVED

%We have to use a negative timespan on the planet approaching hyperbola
%(positive time going backwards); the entire timespan is unknown

%conditions for heliocentric arcs plotting
r_peri = r_peri_modulus.*[0;1;0];
yo = [r_earth+r_peri V_neg]; %arrival state vector
y1 = [r_earth+r_peri V_pos]; %departing state vector
T = 10000000; 
tspano = linspace(20000000,0,10000); %time span where the state vector is defined
tspan1 = linspace(0,20000000,10000);
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_sun), tspan1, y1, options);
[T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_sun), tspano, yo, options);
Y1 = Y1./AU;
Yo = Yo./AU;
figure(3)
grid on
hold on
axis equal
plot3(Yo(:,1),Yo(:,2),Yo(:,3),'b-');
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(0,0,0,'yo');
%% Exercise 1b

%Data definition
mu_planet = astroConstants(13); %Earth case
mu_sun = astroConstants(4);
AU = astroConstants(2);
v_inf_neg = [15.1;0;0]; %km/s
r_earth = [1;0;0]*AU; %Earth's position, km
R_earth = 6371; %Earth's radius, km

%different impact parameters
impact_par_modulus = [1:1:40]'*R_earth;
G = length(impact_par_modulus); %just a dummy variable for sizing

v_inf = norm(v_inf_neg); %modulus of v_inf_neg, km/s
a = -mu_planet/(v_inf^2); %semimajor axis, km

%assumption: u = [0;0;1]
u = [0;0;1];

turning_angle = [];
r_peri_modulus = [];
e = [];
deltaV_modulus = [];
v_inf_pos = [];
V_pos = [];

v_planet_modulus = sqrt(mu_sun/norm(r_earth));
v_planet = v_planet_modulus.*[0;1;0]; %oriented in the Y direction

%computing capital V- (vectors in the heliocentric frame)
V_neg = v_planet+v_inf_neg;

for i=1:G
    turning_angle(i) = 2*atan(-a/impact_par_modulus(i)); %Rad
    e(i) = 1/sin(turning_angle(i)/2); %eccentricity
    deltaV_modulus(i) = 2*v_inf*sin(turning_angle(i)/2); %km/s
    r_peri_modulus(i) = a*(1-e(i)); %pericenter radius, km
    %rotating the approach velocity vector
    v_inf_pos(i,:) = rodrigues(v_inf_neg,u,-turning_angle(i)); %negative angle, leading
    %computing capital V's (vectors in the heliocentric frame)
    V_pos(i,:) = v_planet'+v_inf_pos(i,:);
    
    %conditions for incoming heliocentric arc plotting
    r_peri = r_peri_modulus(i).*[0;1;0];
    yo = [r_earth+r_peri V_neg]; %arrival state vector
    y1 = [r_earth+r_peri V_pos(i,:)']; %departing state vector
    % T = 10000000; 
    tspano = linspace(20000000,0,10000); %time span where the state vector is defined
    tspan1 = linspace(0,20000000,10000);
    %[dy] = twobodyproblem_ode(tspan, y1,muEarth)
    options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
    [T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_sun), tspan1, y1, options);
    [T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_sun), tspano, yo, options);
    Y1 = Y1./AU;
    Yo = Yo./AU;
    figure(4)   
    grid on
    hold on;
    plot3(Y1(:,1),Y1(:,2),Y1(:,3));
end

%adding the Sun and the incoming arc
figure(4)
hold on
axis equal
legend('delta1','delta2','delta3');
plot3(Yo(:,1),Yo(:,2),Yo(:,3));
plot3(0,0,0,'yo');

%%
%Impact parameter vs min flyby altitude
figure(5)
plot(impact_par_modulus./R_earth,(r_peri_modulus')./R_earth)
hold on
plot(impact_par_modulus./R_earth,(turning_angle')*180/pi)
grid on
