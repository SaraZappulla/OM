%% Earth to Ganymed Transfer with saturn flyby

clc 
clear all

%Warning: do not exceed 7 full months on both the departure and arrival
%windows to avoid singularities (null transfer deltaV)
dep_time1 = date2mjd2000([2009,06,18,05,32,0]); %18 june 2009
arr_time1 = date2mjd2000([2009,06,22,12,0,0]); %22 june 2009
mu1 = astroConstants(13); %Earth
mu2 = astroConstants(16); %Saturn

% time_line = [];
% deltavtot = [];
% for i = 1:time_length1
%     for j = 1:time_length2
%         t1 = departing_window(i)*3600*24;
%         t2 = arrival_window(j)*3600*24;
%         [kep1,ksun] = uplanet(departing_window(i), 3); %number 3 corresponds to Earth
%         [kep2,ksun] = uplanet(arrival_window(j), 6); %number 6 corresponds to Saturn
%         [kep,mass,M,d] = ephNEO(arrival_window(j),65);
%         kep1(6) = kep1(6)*180/pi;
%         kep2(6) = kep2(6)*180/pi;
%         kep1(3) = kep1(3)*180/pi;
%         kep2(3) = kep2(3)*180/pi;
%         kep1(4) = kep1(4)*180/pi;
%         kep2(4) = kep2(4)*180/pi;
%         kep1(5) = kep1(5)*180/pi;
%         kep2(5) = kep2(5)*180/pi;
%         deltavtot(i,j) = deltavtot_computation(t1,t2,kep1,kep2,mu1,mu2,ksun);
%         time_line(i,j) = arrival_window(j) - departing_window(i);
% 
%         
%     end
%    
% end

%% optimisation graph

% to obtain dates inside the graph
dep_wind = [];
arr_wind = [];
%for i = 1:time_length1
%        dep_wind(i,:)= mjd20002date(departing_window(i))
%        dep_windstr(i,1) = num2str(dep_wind(i,1));
%end
%for j = 1:time_length2
%        arr_wind(j,:)= mjd20002date(arrival_window(j))
%        arr_windstr(j,1) = num2str(arr_wind(j,1));
%end
mindelta = min(min(deltavtot));
maxdelta = max(max(deltavtot));

figure(1)
hold on
contour(departing_window,arrival_window,deltavtot',linspace(mindelta-1,maxdelta+1,100));
c = colorbar;
c.Label.String = '\DeltaV[km/s]';
%z = contour(departing_window,arrival_window,time_line',[60,120,180,240,300],'k');
%d = colorbar(z,'off')
xlabel('Departing Date')
ylabel('Arrival Date')

%%


for i=1:time_length1
    for j = 1:time_length2
        if deltavtot(i,j) == mindelta
            k=i
            c=j
        end
    end
end

cheap_departure = departing_window(k);
cheap_intermidiate = arrival_window(c);
cheap_departure_date = mjd20002date(cheap_departure);
cheap_intermidiate_date = mjd20002date(cheap_intermidiate);

        [kep1_min,muSun] = uplanet(cheap_departure, 3);
        [kep2_min,muSun] = uplanet(cheap_intermidiate, 6);

        kep1_min(6) = kep1_min(6)*180/pi;
        kep2_min(6) = kep2_min(6)*180/pi;
                kep1_min(3) = kep1_min(3)*180/pi;
        kep2_min(3) = kep2_min(3)*180/pi;
                kep1_min(4) = kep1_min(4)*180/pi;
        kep2_min(4) = kep2_min(4)*180/pi;
                kep1_min(5) = kep1_min(5)*180/pi;
        kep2_min(5) = kep2_min(5)*180/pi;

%computing all the orbits
%orbit 1
T1 = 2*pi * sqrt(kep1_min(1)^3/muSun);
tspan1 = linspace(0,T1,10000);
[r1,v1]=kepl_to_car(kep1_min(1),kep1_min(2),kep1_min(3),kep1_min(4),kep1_min(5),kep1_min(6),muSun);
y1 = [r1 v1];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t1,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan1, y1, options);

%orbit 2
T2 = 2*pi * sqrt(kep2_min(1)^3/muSun);
tspan2 = linspace(0,T2,10000);
[r2,v2]=kepl_to_car(kep2_min(1),kep2_min(2),kep2_min(3),kep2_min(4),kep2_min(5),kep2_min(6),muSun);
y2 = [r2 v2];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t2,Y2] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan2, y2, options);

%transfer orbit
ToF = (cheap_intermidiate - cheap_departure)*24*3600;
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muSun, 0, 0, 0 );
TT = 2*pi * sqrt(A^3/muSun);
tspanT = linspace(0,TT,10000);
yT = [r1' VI];
options= odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[tT,YT] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanT, yT, options);
tspanTOF = linspace(0,ToF,10000);
[tT,Yarc] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanTOF, yT, options);

%%
figure(1)
hold on
plot3(departing_window(k),arrival_window(c),deltavtot(k,c),'o')

figure(2)
hold on
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');

plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-');
plot3(Yarc(:,1),Yarc(:,2),Yarc(:,3),'c-');
plot3(YT(:,1),YT(:,2),YT(:,3),'c--');
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'bo');
legend('initial orbit','final orbit','transfer arc','transfer orbit','initial transfer point','final transfer point','Location','northeastoutside')



%% section 2: Heliocentric orbit Earth to Saturn
k = 16; %choosen depending on minimum speed
mintime = k*365; %minimum transfer possible time to get to sarutn
%intermidiate_window represent all teh possible time to arrive to saturn
%from the starting from earth


deltav1 = [];
deltav2 = [];
deltavtot = [];
for i = 1: time_length1
    intermidiate_window = linspace(departing_window(i)+mintime,departing_window(i)+mintime+time_length3,time_length3);
    for j = 1: time_length3
        % computing planets ephemerides
        [kep1,muSun] = uplanet(departing_window(i), 3); %number 3 corresponds to Earth
        [kep2,muSun] = uplanet(intermidiate_window(j), 6); %number 6 corresponds to Saturn
        kep1(3:6) = kep1(3:6)*180/pi;
        kep2(3:6) = kep2(3:6)*180/pi;
        
        % computing delta v of the heliocentric leg
        t1 = departing_window(i)*3600*24; %seconds
        t2 = intermidiate_window(j)*3600*24; %seconds
        ToF = t2-t1; % time that the explorer takes to fly on the lambert arc [s]
        [r1,v1]=kepl_to_car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),muSun);
        [r2,v2]=kepl_to_car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muSun);
        [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1, r2, ToF, muSun, 0, 0, 0 );
        kep_t = [A,E,0,0,0,THETA]; %transfer orbit parameters
        
        % VI and VF s/c velocity in heliocentric
        % v1 and v2 of the two planets
        deltav1(i,j) = norm(VI'-v1); %v infinite plus of the spacecraft in Earth SOI
        deltav2(i,j) = norm(v2-VF'); %v infinite minus of the spacecraft in Saturn SOI
        deltavtot(i,j) = (deltav1(i,j))+(deltav2(i,j));
    end
end

minv = min(min(deltavtot))
% heliocentric leg plot of minimum delta v

for i = 1: time_length1
    intermidiate_window = linspace(departing_window(i)+mintime,departing_window(i)+mintime+time_length3,time_length3);
    for j = 1: time_length3
        if deltavtot(i,j) == minv
            k=i
            c=j
            inter_window = intermidiate_window;
        end
    end
end

cheap_departure = departing_window(k);
cheap_intermidiate = inter_window(c);
cheap_departure_date = mjd20002date(cheap_departure);
cheap_intermidiate_date = mjd20002date(cheap_intermidiate);

% computation of algebraic ephemerides of the planets
[kep1_min,muSun] = uplanet(cheap_departure, 3);
[kep2_min,muSun] = uplanet(cheap_intermidiate, 6);

% transformation from radiants to degrees
kep1_min(6) = kep1_min(6)*180/pi;
kep2_min(6) = kep2_min(6)*180/pi;
kep1_min(3) = kep1_min(3)*180/pi;
kep2_min(3) = kep2_min(3)*180/pi;
kep1_min(4) = kep1_min(4)*180/pi;
kep2_min(4) = kep2_min(4)*180/pi;
kep1_min(5) = kep1_min(5)*180/pi;
kep2_min(5) = kep2_min(5)*180/pi;

% propagation of all the orbits
%orbit 1
T1 = 2*pi * sqrt(kep1_min(1)^3/muSun);
tspan1 = linspace(0,T1,10000);
[r1,v1]=kepl_to_car(kep1_min(1),kep1_min(2),kep1_min(3),kep1_min(4),kep1_min(5),kep1_min(6),muSun);
y1 = [r1 v1];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t1,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan1, y1, options);

%orbit 2
T2 = 2*pi * sqrt(kep2_min(1)^3/muSun);
tspan2 = linspace(0,T2,10000);
[r2,v2]=kepl_to_car(kep2_min(1),kep2_min(2),0,kep2_min(4),kep2_min(5),kep2_min(6),muSun);
y2 = [r2 v2];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t2,Y2] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan2, y2, options);
v_planet = v2;

%transfer orbit
ToF = (cheap_intermidiate - cheap_departure)*24*3600;
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muSun, 0, 0, 0 );
TT = 2*pi * sqrt(A^3/muSun);
tspanT = linspace(0,TT,10000);
yT = [r1' VI];
options= odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[tT,YT] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanT, yT, options);
tspanTOF = linspace(0,ToF,10000);
[tT,Yarc] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanTOF, yT, options);
v_inf_neg_s = v2-VF';
r_saturn = r2;
V_neg_s = VF;
% heliocentric orbits plot

figure(1)
hold on
grid on
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-');
plot3(Yarc(:,1),Yarc(:,2),Yarc(:,3),'c-');
plot3(YT(:,1),YT(:,2),YT(:,3),'c--');
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'bo');
legend('initial orbit','final orbit','transfer arc','transfer orbit','initial transfer point','final transfer point','Location','northeastoutside')



%% section 1: Earth sphere of influence

% circular orbit
% put a,e,...
R_earth = 6378.1;
a = 6778.1;
e = 0;
vcircular = sqrt(mu1/a);

% hyperbolic orbit
rp = a; %pericener of parabola is on the circular orbit
v_infinite_plus = deltav1(k,c);
% conservation of energy to get v at pericenter
vp = sqrt(v_infinite_plus^2 + 2*mu1/rp);

%characterisation of hyperbola
impact_par_modulus = vp * rp / v_infinite_plus;
ahyp = -mu1/(v_infinite_plus^2); %semimajor axis, km
turning_angle = 2*atan(-ahyp/impact_par_modulus); %Rad
ehyp = 1/sin(turning_angle/2); %eccentricity
%r_peri_modulus = a*(1-e); %pericenter radius, km
alpha = 2*pi-turning_angle/2;

i = 0;
OM = 0;
om = 0;
theta = 0;

deltav_manoeuvre = abs(vp-vcircular);

%% orbit propagation
%circular orbit
r_peri_c = rp.*[1;0;0]; 
r_peri_hyp = rp.*[cos(alpha);-sin(alpha);0];
yo = [r_peri_c vcircular.*[0;1;0]]; %circular orbit state vector (pericenter)
y1 = [r_peri_hyp vp.*[sin(alpha);cos(alpha);0]]; %departing hyperbola state vector (pericenter)
%T = 10000000; 
tspano = linspace(10000,0,10000); %time span where the state vector is defined
tspan1 = linspace(0,1000,10000);
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu1), tspan1, y1, options);
[T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu1), tspano, yo, options);
%Yo = Yo./R_earth;
%Y1 = Y1./R_earth;
figure(2);
grid on
hold on
axis equal
earth_sphere
plot3(Yo(:,1),Yo(:,2),Yo(:,3),'b-');
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');

%% section 4: Saturn to Ganymed

t = 12; %choosen depending on minimum speed
mintime = t*365; %minimum transfer possible time to get to sarutn
%intermidiate_window represent all teh possible time to arrive to saturn
%from the starting from earth

deltav1_s = [];
deltav2_s = [];
deltavtot_s = [];
for i = 1: time_length2
    intermidiate_window = linspace(arrival_window(end)-mintime-time_length3,arrival_window(end)-mintime,time_length3);
    for j = 1: time_length3
        % computing planets ephemerides
        [kep1,mass,M] = ephNEO(arrival_window(i),65);
        [kep2,muSun] = uplanet(intermidiate_window(j), 6); %number 6 corresponds to Saturn
        kep1(6) = kep1(6)*180/pi;
        kep2(6) = kep2(6)*180/pi;
        kep1(3) = kep1(3)*180/pi;
        kep2(3) = kep2(3)*180/pi;
        kep1(4) = kep1(4)*180/pi;
        kep2(4) = kep2(4)*180/pi;
        kep1(5) = kep1(5)*180/pi;
        kep2(5) = kep2(5)*180/pi;
        
        % computing delta v of the heliocentric leg
        t1 = arrival_window(i)*3600*24; %seconds
        t2 = intermidiate_window(j)*3600*24; %seconds
        ToF = t2-t1; % time that the explorer takes to fly on the lambert arc [s]
        [r1,v1]=kepl_to_car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),muSun);
        [r2,v2]=kepl_to_car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muSun);
        [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1, r2, ToF, muSun, 0, 0, 0 );
        kep_t = [A,E,0,0,0,THETA]; %transfer orbit parameters
        
        % VI and VF s/c velocity in heliocentric
        % v1 and v2 of the two planets
        deltav1_s(i,j) = norm(VI'-v1); %v infinite plus of the spacecraft in Earth SOI
        deltav2_s(i,j) = norm(v2-VF'); %v infinite minus of the spacecraft in Saturn SOI
        deltavtot_s(i,j) = (deltav1(i,j))+(deltav2(i,j));
    end
end

minv_s = min(min(deltavtot_s))
% heliocentric leg plot of minimum delta v

for i = 1: time_length2
    intermidiate_window = linspace(arrival_window(end)-mintime-time_length3,arrival_window(end)-mintime,time_length3);
    for j = 1: time_length3
        if deltavtot_s(i,j) == minv_s
            k=i
            c=j
            inter_window = intermidiate_window;
        end
    end
end

cheap_arrival = arrival_window(k);
cheap_intermidiate = inter_window(c);
cheap_arrival_date = mjd20002date(cheap_arrival);
cheap_intermidiate_date = mjd20002date(cheap_intermidiate);

% computation of algebraic ephemerides of the planets
[kep1_min,mass,M] = ephNEO(cheap_arrival,65);
[kep2_min,muSun] = uplanet(cheap_intermidiate, 6);

% transformation from radiants to degrees
kep1_min(6) = kep1_min(6)*180/pi;
kep2_min(6) = kep2_min(6)*180/pi;
kep1_min(3) = kep1_min(3)*180/pi;
kep2_min(3) = kep2_min(3)*180/pi;
kep1_min(4) = kep1_min(4)*180/pi;
kep2_min(4) = kep2_min(4)*180/pi;
kep1_min(5) = kep1_min(5)*180/pi;
kep2_min(5) = kep2_min(5)*180/pi;

% propagation of all the orbits
%orbit 1
T1 = 2*pi * sqrt(kep1_min(1)^3/muSun);
tspan1 = linspace(0,T1,10000);
[r1,v1]=kepl_to_car(kep1_min(1),kep1_min(2),kep1_min(3),kep1_min(4),kep1_min(5),kep1_min(6),muSun);
y1 = [r1 v1];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t1,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan1, y1, options);

%orbit 2
T2 = 2*pi * sqrt(kep2_min(1)^3/muSun);
tspan2 = linspace(0,T2,10000);
[r2,v2]=kepl_to_car(kep2_min(1),kep2_min(2),0,kep2_min(4),kep2_min(5),kep2_min(6),muSun);
y2 = [r2 v2];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t2,Y2] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan2, y2, options);

%transfer orbit
ToF = (cheap_arrival - cheap_intermidiate)*24*3600;
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muSun, 0, 0, 0 );
TT = 2*pi * sqrt(A^3/muSun);
tspanT = linspace(0,TT,10000);
yT = [r1' VI];
options= odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[tT,YT] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanT, yT, options);
tspanTOF = linspace(0,ToF,10000);
[tT,Yarc] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanTOF, yT, options);
V_pos_s = VI;
% heliocentric orbits plot

figure(1)
hold on
grid on
axis equal
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-');
plot3(Yarc(:,1),Yarc(:,2),Yarc(:,3),'c-');
plot3(YT(:,1),YT(:,2),YT(:,3),'c--');
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'bo');
legend('initial orbit','final orbit','transfer arc','transfer orbit','initial transfer point','final transfer point','Location','northeastoutside')



%% section 3: Saturn flyby

%Data definition
mu_planet = astroConstants(16); %Saturn case
mu_sun = astroConstants(4);
AU = astroConstants(2);
%impact_par_modulus = 9200; %km
%r_earth = [1;0;0]*AU; %Earth's position, km
%v_inf_neg = deltav2(k,c); %take from the first heliocentric leg
R_saturn = astroConstants(26);
R_rings = 135000; %km

% definition of impact parameter
impact_par_modulus = [1.01:0.01:1.2]*R_saturn; %probably we need a smaller one like 8000
G = length(impact_par_modulus); %just a dummy variable for sizing
%v_inf_neg = v_inf_neg*u;

v_inf_neg = norm(v_inf_neg_s); %modulus of v_inf_neg_s, km/s
a = -mu_planet/(v_inf_neg^2); %semimajor axis, km

turning_angle = [];
r_peri_modulus = [];
e = [];
deltaV_modulus = [];
v_inf_pos = [];
V_pos = [];

v_planet_modulus = norm(v_planet);
%computing capital V- (vectors in the heliocentric frame)
V_neg = v_planet+v_inf_neg_s; 

%u is the unit vector around which the velocity turn of an angle delta. It
%is considered here as the normal to the lambert plane (in heliocentric
%system) 
u = cross(r_saturn,V_neg)/norm(cross(r_saturn,V_neg));

for i=1:G
    turning_angle(i) = 2*atan(-a/impact_par_modulus(i)); %Rad
    e(i) = 1/sin(turning_angle(i)/2); %eccentricity
    deltaV_modulus(i) = 2*v_inf_neg*sin(turning_angle(i)/2); %km/s
    r_peri_modulus(i) = a*(1-e(i)); %pericenter radius, km
    %rotating the approach velocity vector
    v_inf_pos(i,:) = rodrigues(v_inf_neg_s,u,-turning_angle(i)); %negative angle, leading
    %computing capital V's (vectors in the heliocentric frame)
    V_pos(i,:) = v_planet'+v_inf_pos(i,:);
    beta = (pi-turning_angle)/2;
    
    v_neg_direction = v_inf_neg_s/(v_inf_neg);
    rp_direction = rodrigues(v_inf_neg_s,v_neg_direction,-beta);
    rp_direction_unit = rp_direction/(norm(rp_direction));
    %conditions for incoming heliocentric arc plotting
    r_peri = r_peri_modulus(i).*rp_direction_unit;
    yo = [r_peri; v_inf_neg_s]; %arrival state vector
    y1 = [r_peri; v_inf_pos(i,:)']; %departing state vector
    % T = 10000000; 
    tspano = linspace(200000,0,10000); %time span where the state vector is defined
    tspan1 = linspace(0,200000,10000);
    %[dy] = twobodyproblem_ode(tspan, y1,muEarth)
    options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
    [T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_planet), tspan1, y1, options);
    [T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_planet), tspano, yo, options);
    Y1 = Y1./AU;
    Yo = Yo./AU;
    figure(6)   
    grid on
    hold on;
    plot3(Y1(:,1),Y1(:,2),Y1(:,3));
    i
end

%% Flyby not natural

%V_pos_s from the Ganymed side
%V_neg_s from the Eart side
%Data definition
mu_planet = astroConstants(16); %Saturn case
mu_sun = astroConstants(4);
AU = astroConstants(2);
G = astroConstants(1);
%r_earth = [0;-1;0]*AU; %Earth's position, km
%R_earth = 6371; %Earth's radius, km
%h_atm = 100; %Earth's atmosphere altitude, km 
%computing Vplanet
r_saturn = astroConstants(26); % in planetocentric
r_rings = 135000; %km in planetocentric

%calculating approach and departure velocities
v_inf_neg = V_neg_s'-v_planet;
v_inf_pos = V_pos_s'-v_planet;

%calculating the turning angle
v_neg = norm(v_inf_neg);
v_pos = norm(v_inf_pos);

turning_angle = acos(dot(v_inf_neg,v_inf_pos)/(v_pos*v_neg));

%calculating r_peri
turning_angle_function = @(r_peri) (2*asin(1/(1+r_peri*(v_neg^2)/mu_planet)))/2+...
    (2*asin(1/(1+r_peri*(v_pos^2)/mu_planet)))/2;%-turning_angle;

%definition of an initial guess
r_peri_guess = r_saturn+300; %Saturn's radius + random height, km
r_peri_fzero = fsolve(turning_angle_function,r_peri_guess);
if r_peri_fzero < (r_saturn+200)
    fprintf('error, too low pericenter radius\n');
end
m_saturn = mu_planet/G;
m_sun = mu_sun/G;
R_saturn_mean = 9.5826*AU;
SOI = R_saturn_mean * (m_saturn/m_sun)^0.4;

if r_peri_fzero > (SOI)
    fprintf('error, too high pericenter radius\n');
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
deltaV_flyby = norm(V_pos_s-V_neg_s);

%arcs plotting, planetocentric
r_peri = r_peri_fzero.*[1;0;0]; 
yo = [r_peri v_peri_in.*[0;1;0]]; %arrival state vector (pericenter)
y1 = [r_peri v_peri_out.*[0;1;0]]; %departing state vector (pericenter)
%T = 10000000; 
tspano = linspace(1000000,0,100000); %time span where the state vector is defined
tspan1 = linspace(0,1000000,100000);
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_planet), tspan1, y1, options);
[T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_planet), tspano, yo, options);
figure(5);
grid on
hold on
axis equal
plot3(Yo(:,1),Yo(:,2),Yo(:,3),'b-');
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(0,0,0,'.k','MarkerSize',50)
% impact parameter is a DOF maybe


%% Total problem

% 1. understant which is the synodic period of the 3 bodies and how often
% the combinations get repeated (useful to understand which should be the
% pace considered)
% 2. Compute the 2 legs in any time combination considering:
%       - time of departure 
%       - time of flyby
%       - time of arrival
% 3. Do porkchop plot to undertand which manoeuvre is the less expensive

%% computing synodic period
k = 16;
mu_Sun = astroConstants(4);
[kep1,muSun] = uplanet(departing_window(1), 3); %number 3 corresponds to Earth
[kep2,muSun] = uplanet(arrival_window(1)+k*365, 6); %number 6 corresponds to Saturn
[kep3,mass,M] = ephNEO(arrival_window(1),65); % M is the mean anomaly of asteroid

mu1 = astroConstants(13);
T_earth = 2*pi* sqrt((kep1(1))^3/muSun); % mean motion of the earth
n1 = 2*pi / (T_earth)

mu2 = astroConstants(16);
T_saturn = 2*pi* sqrt((kep2(1))^3/muSun); % mean motion of the earth
n2 = 2*pi / (T_saturn)

mu3 = G*mass;
T_ganymed = 2*pi* sqrt((kep3(1))^3/muSun); % mean motion of the earth
n3 = 2*pi / (T_ganymed)

syn_se = 2*pi / abs(n2-n1); %earth-saturn synodic period
syn_se_days = syn_se /3600/24/365

syn_as = 2*pi / abs(n3-n2); %saturn_ganymed synodic period
syn_as_days = syn_as /3600/24/365

syn_ae = 2*pi / abs(n3-n1); %saturn_ganymed synodic period
syn_ae_days = syn_ae /3600/24/365

total_syn = (syn_se_days * syn_as_days * syn_ae_days)

%% calculating the min delta v
% years from:
deltav1 = [];
deltav2 = [];
deltavtot_es = [];
deltav1_s = [];
deltav2_s = [];
deltavtot_sg = [];
deltayears = [];
k = 1;
c=1;
alpha = [4000:100:9125]; % earth to saturn minimum time
gamma = [1825:100:7300]; %saturn to ganymed

for alphai = 1:length(alpha) %assuming you are taking 10 to 25 years to arrive to saturn 
    deltavtot_es = 100;
    for i = 1:time_length1 %considering time window to leave earth 5 years span
        intermidiate_window = linspace(departing_window(i)+alpha(alphai),departing_window(i)+alpha(alphai)+time_length3,time_length3);
        deltavpreci = deltavtot_es
        for j = 1:time_length3 % consdiering time window to arrive to saturn
            % put earth-saturn leg
            [kep1,muSun] = uplanet(departing_window(i), 3); %number 3 corresponds to Earth
            [kep2,muSun] = uplanet(intermidiate_window(j), 6); %number 6 corresponds to Saturn
            kep1(3:6) = kep1(3:6).*(180/pi);
            kep2(3:6) = kep2(3:6).*(180/pi);
            
            % computing delta v of the heliocentric leg
            t1 = departing_window(i)*3600*24; %seconds
            t2 = intermidiate_window(j)*3600*24; %seconds
            ToF = t2-t1; % time that the explorer takes to fly on the lambert arc [s]
            [r1,v1]=kepl_to_car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),muSun);
            [r2,v2]=kepl_to_car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muSun);
            [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1, r2, ToF, muSun, 0, 0, 0 );
            kep_t = [A,E,0,0,0,THETA]; %transfer orbit parameters
            
            deltavprecj = deltavtot_es;
            % VI and VF s/c velocity in heliocentric
            % v1 and v2 of the two planets
%             deltav1(i,j,c+1,k+1) = norm(VI'-v1); %v infinite plus of the spacecraft in Earth SOI
%             deltav2(i,j,c+1,k+1) = norm(v2-VF'); %v infinite minus of the spacecraft in Saturn SOI
%             deltavtot_es(i,j,c,k) = (deltav1(i,j,c+1,k+1))+(deltav2(i,j,c+1,k+1));

            deltav1 = norm(VI'-v1); %v infinite plus of the spacecraft in Earth SOI
            deltav2= norm(v2-VF'); %v infinite minus of the spacecraft in Saturn SOI
            deltavtot_es = (deltav1)+(deltav2);
            deltavnowj = deltavtot_es;

            if deltavprecj > deltavnowj
                deltayears(alphai) = deltavnowj;
            else
                deltayears(alphai) = deltavprecj;
            end

            %earth side already in the cost?
            
    for gammai = 1:length(gamma) %assuming you are taking 5 to 20 years to arrive to ganymed 
        deltavtot_sg = 100;
        for k = 1:time_length2 % considering time window to arrive to ganymed
            intermidiate_window_2 = linspace(arrival_window(end)-gamma(gammai)-time_length3,arrival_window(end)-gamma(gammai),time_length3);
            deltavpreck = deltavtot_sg;
            for c = 1:time_length3 % considering time window to leave saturn or time of flyby
                %put saturn-ganyned leg

                [kep3,mass,M] = ephNEO(arrival_window(k),65);
                [kep2,muSun] = uplanet(intermidiate_window_2(c), 6); %number 6 corresponds to Saturn
                kep3(3:6) = kep3(3:6).*(180/pi);
                kep2(3:6) = kep2(3:6).*(180/pi);
                
                % computing delta v of the heliocentric leg
                t3 = arrival_window(k)*3600*24; %seconds
                t2 = intermidiate_window_2(c)*3600*24; %seconds
                ToF = t3-t2; % time that the explorer takes to fly on the lambert arc [s]
                [r3,v3]=kepl_to_car(kep3(1),kep3(2),kep3(3),kep3(4),kep3(5),kep3(6),muSun);
                [r2,v2]=kepl_to_car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muSun);
                [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r2, r3, ToF, muSun, 0, 0, 0 );
                kep_t = [A,E,0,0,0,THETA]; %transfer orbit parameters
                
                deltavprecc = deltavtot_sg;
                % VI and VF s/c velocity in heliocentric
                % v1 and v2 of the two planets
                %deltav1_s(i,j,c,k) = norm(VI'-v1); %v infinite plus of the spacecraft in Earth SOI
                %deltav2_s(i,j,c,k) = norm(v2-VF'); %v infinite minus of the spacecraft in Saturn SOI
                %deltavtot_sg(i,j,c,k) = (deltav1_s(i,j,c,k))+(deltav2_s(i,j,c,k));
                deltav1 = norm(VI'-v1); %v infinite plus of the spacecraft in Earth SOI
                deltav2= norm(v2-VF'); %v infinite minus of the spacecraft in Saturn SOI
                deltavtot_sg = (deltav1)+(deltav2);
                deltavnowc = deltavtot_sg;

                if deltavprecc > deltavnowc
                    deltayears(alphai,gammai) = deltavnowc;
                else
                    deltayears(alphai,gammai) = deltavprecc;
                end
                %put fly by

            end
        deltavnowk = deltavtot_sg; 
        if deltavpreck > deltavnowk
                    deltayears(alphai,gammai) = deltavnowk;
                else
                    deltayears(alphai,gammai) = deltavpreck;
                end
        end
        
    end
    k
        end
        deltavnowi = deltavtot_es;
        if deltavprec > deltavnow
            deltayears(alphai) = deltavnowi;
        else
            deltayears(alphai) = deltavpreci;
        end

    end
    j
end

%deltavtot = deltav_es + deltav_sg;
figure(1)
hold on
contour(alpha,beta,deltavyears);%,linspace(mindelta-1,maxdelta+1,100)
c = colorbar;
c.Label.String = '\DeltaV[km/s]';
%z = contour(departing_window,arrival_window,time_line',[60,120,180,240,300],'k');
%d = colorbar(z,'off')
xlabel('Departing Date')
ylabel('Arrival Date')

                

