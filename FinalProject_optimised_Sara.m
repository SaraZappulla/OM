%% Total problem

%OUR ASTEROID IS 1994GK, NOT GANYMED

% 1. problem data
% 2. understant which is the synodic period of the 3 bodies and how often
% the combinations get repeated (useful to understand which should be the
% time window considered)
% 3. Compute the 2 legs in any time combination considering:
%       - time of departure 
%       - time of flyby
%       - time of arrival
% 4. Do porkchop plot to undertand which manoeuvre is the less expensive
% 5. propagate cheaper manoeuvres
% 6. compute earth parking orbit
% 7. compute saturn fly by


%% 1. problem data


dep_time1 = date2mjd2000([2030,11,30,12,0,0]); %30 november 2030
dep_time2 = date2mjd2000([2035,11,29,12,0,0]); %29 november 2035
arr_time1 = date2mjd2000([2060,10,29,12,0,0]); %29 october 2060
arr_time2 = date2mjd2000([2065,10,28,12,0,0]); %28 october 2065
time_length1 = 60; %it gives a pace of a month but is an entire number
time_length_days = (dep_time2 - dep_time1); 
time_length2 = 60; %(arr_time2 - arr_time1)/30; it gives a pace of a month but is an entire number
time_length3 = 60;
departing_window = linspace(dep_time1,dep_time2,time_length1);
arrival_window = linspace(arr_time1,arr_time2,time_length2);
mu1 = astroConstants(13); %Earth
mu2 = astroConstants(16); %Saturn

%% 2.computing synodic period
k = 16; % random number to compute the flyby time it correspond to years

% bodies ephemerides
[kep1,muSun] = uplanet(departing_window(1), 3); %number 3 corresponds to Earth
[kep2,muSun] = uplanet(arrival_window(1)+k*365, 6); %number 6 corresponds to Saturn
[kep3,mass,M] = ephNEO(arrival_window(1),65); % M is the mean anomaly of asteroid
G = astroConstants(1); %universal gravitational constant

%earth mean motion
mu1 = astroConstants(13);
T_earth = 2*pi* sqrt((kep1(1))^3/muSun); 
n1 = 2*pi / (T_earth);

%saturn mean motion
mu2 = astroConstants(16);
T_saturn = 2*pi* sqrt((kep2(1))^3/muSun);
n2 = 2*pi / (T_saturn);

%ganymed mean motion
mu3 = G*mass;
T_ganymed = 2*pi* sqrt((kep3(1))^3/muSun); 
n3 = 2*pi / (T_ganymed);

syn_se = 2*pi / abs(n2-n1); %earth-saturn synodic period
syn_se_years = syn_se /3600/24/365

syn_as = 2*pi / abs(n3-n2); %saturn_ganymed synodic period
syn_as_years = syn_as /3600/24/365

syn_ae = 2*pi / abs(n3-n1); %saturn_ganymed synodic period
syn_ae_years = syn_ae /3600/24/365

% total synodic period of the three bodies in years
total_syn = (syn_se_years * syn_as_years * syn_ae_years)



%% 3. calculating the delta v for every combination

deltavtot = []; %used for total optimisation
deltayears = []; %used for pork chop plot

%alpha must be higher than 1800 days since it's not possible to do the
%maneuver due to requirements in velocity too high - this can be verified
%by finding the minimum transfer time on ToF_i_alphai [s]

%Earth to Saturn minimum time (days) = 1800 - computed with ToF_i_alphai 
%Saturn to Ganymed maximum time (days) = 9126 - computed with ToF_k_alphai

alpha = [1800:100:7000]; % earth to saturn minimum time 
%ToF_i_alphai = [];
%ToF_k_alphai = []; %to check that the final element of the intermediate window allows the spacecraft to reach Ganymed before the last arrival date
    
for i = 1:time_length1 %considering time window to leave earth 5 years span
    for k = 1:time_length2 % considering time window to arrive to ganymed
        deltavprec = 1000;
        for alphai = 1:length(alpha) %assuming you are taking 10 to 25 years to arrive to saturn
            % earth-saturn leg
            intermediate_window(alphai) = (departing_window(i)+alpha(alphai));
            [kep1,muSun] = uplanet(departing_window(i), 3); %number 3 corresponds to Earth
            [kep2,muSun] = uplanet(intermediate_window(alphai), 6); %number 6 corresponds to Saturn
            [kep3,mass,M] = ephNEO(arrival_window(k),65); %Ganymed ephemerides

            kep1(3:6) = kep1(3:6).*(180/pi);
            kep2(3:6) = kep2(3:6).*(180/pi);
            kep3(3:6) = kep3(3:6).*(180/pi);
           
            % computing delta v of the heliocentric leg 1
            t1 = departing_window(i)*3600*24; %seconds
            t2 = intermediate_window(alphai)*3600*24; %seconds
            ToF1 = t2-t1; % time that the explorer takes to fly on the lambert arc [s]
            % ToF_i_alphai(i,alphai) = t2-t1; 
            [r1,v1]=kepl_to_car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),muSun);
            [r2,v2]=kepl_to_car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muSun);
            [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1, r2, ToF1, muSun, 0, 0, 0 );
           
            % VI and VF s/c velocity in heliocentric
            % v1 and v2 of the two planets

            deltav1 = norm(VI-v1'); %v infinite plus of the spacecraft in Earth SOI
            deltav2= norm(v2'-VF); %v infinite minus of the spacecraft in Saturn SOI
            deltavtot_es = (deltav1)+(deltav2);
         
            %saturn-ganyned leg
            
            % computing delta v of the heliocentric leg
            t3 = arrival_window(k)*3600*24; %seconds
            %t2 already defined above
            ToF2 = t3-t2; % time that the explorer takes to fly on the lambert arc [s]
            % ToF_k_alphai(k,alphai) = t3-t2;
            [r3,v3]=kepl_to_car(kep3(1),kep3(2),kep3(3),kep3(4),kep3(5),kep3(6),muSun);
            [r2,v2]=kepl_to_car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muSun);
            [A,P,E,ERROR,VI2,VF2,TPAR,THETA] = lambertMR(r2, r3, ToF2, muSun, 0, 0, 0 );
            % [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r2, r3, ToF_k_alphai(k,alphai), muSun, 0, 0, 0 );
            kep_t = [A,E,0,0,0,THETA]; %transfer orbit parameters

            % VI and VF s/c velocity in heliocentric
            % v1 and v2 of the two planets
            deltav3 = norm(v2-VI2'); %v infinite plus of the spacecraft in Earth SOI
            deltav4= norm(VF2'-v3); %v infinite minus of the spacecraft in Saturn SOI
            deltavtot_sg = (deltav3)+(deltav4);

            deltavnow = deltavtot_sg + deltavtot_es;
            deltavtot(i,k,alphai) = deltavnow;
            if deltavprec > deltavnow
                deltayears(i,k) = deltavnow;
            else
                deltayears(i,k) = deltavprec;
            end
            deltavprec = deltayears(i,k);

            V_neg_s = VF;
            V_pos_s = VI2;

            v_planet = v2; %velocity of saturn
            v_inf_neg = V_neg_s'-v_planet;
            v_inf_pos = V_pos_s'+v_planet;

            %calculating the turning angle
            v_neg = norm(v_inf_neg);
            v_pos = norm(v_inf_pos);

            turning_angle = acos(dot(v_inf_neg,v_inf_pos)/(v_pos*v_neg));

            %calculating r_peri
            turning_angle_function = @(r_peri) +(asin(1/(1+r_peri*(v_neg^2)/mu_planet)))+...
                 (asin(1/(1+r_peri*(v_pos^2)/mu_planet)))-turning_angle

            %definition of an initial guess
            r_peri_guess = r_saturn*10; %Saturn's radius + random height, km
            r_peri_fzero(i,k,alphai) = fzero(turning_angle_function,r_peri_guess);
        end
    end
end

% Unoptimized intermediate window computation
intermediate_window_1 = mjd20002date(departing_window(1)+alpha(1));
intermediate_window_2 = mjd20002date(departing_window(end)+alpha(end));


       
%% 4. porkchop plot represantation

dep_date = [];
arr_date = [];
for i=1:time_length1
    dep_date(i,:) = mjd20002date(departing_window(i));
end
for i=1:time_length2
    arr_date(i,:) = mjd20002date(arrival_window(i));
end

A = num2str(dep_date(:,1:3)); %yymmdd
B = num2str(arr_date(:,1:3));

figure(1)
hold on
contour(departing_window,arrival_window,deltayears',[30:2:70])%linspace(mindelta-1,maxdelta+1,100));
c = colorbar;
c.Label.String = '\DeltaV[km/s]';
xticklabels(A);
yticklabels(B);
xlabel('Departing dates')
ylabel('Arrival dates')
title('Porkchop plot representation')

%% 
t = zeros(time_length1,time_length2,length(alpha))
for i = 1:time_length1 %considering time window to leave earth 5 years span
    for k = 1:time_length2 % considering time window to arrive to ganymed
        deltavprec = 1000;
        for alphai = 1:length(alpha) %assuming you are taking 10 to 25 years to arrive to saturn
            if abs(r_peri_fzero) < r_saturn
                v(i,k,alphai) = 0;
            else
                v(i,k,alphai) = 1;
            end
            if v(i,k,alphai) == 1;
                t(i,k) = 1;
            end
        end

    end
end
%figure(2)
%contour(departing_window,arrival_window,r_peri_fzero')

%% plot surface considering all the 3 data

% figure(2)
% hold on
% plot3(deltavtot(:,1,1),deltavtot(1,:,1),deltavtot(1,1,:)); %linspace(mindelta-1,maxdelta+1,100))
% %c = colorbar;
% %c.Label.String = '\DeltaV[km/s]';
% %z = contour(departing_window,arrival_window,time_line',[60,120,180,240,300],'k');
% %d = colorbar(z,'off')
% xlabel('Duration for heliocentric leg 1 [days]')
% ylabel('Duration for heliocentric leg 2 [days]')

%% 5. optimised transfer 
mindelta = min(min(min(deltavtot)))

for i=1:time_length1
    for j = 1:time_length2
        for alphai = 1:length(alpha)
            if deltavtot(i,j,alphai) == mindelta
                k=i
                d=j
                gamma = alphai
            end
        end
    end
end

cheap_departure = departing_window(k);
cheap_arrival = arrival_window(d);
cheap_intermediate = departing_window(k) + alpha(gamma);
cheap_departure_date = mjd20002date(cheap_departure);
cheap_arrival_date = mjd20002date(cheap_arrival);
cheap_intermediate_date = mjd20002date(cheap_intermediate);

%% 5. heliocentric leg 1

% computation of algebraic ephemerides of the planets
[kep1_min,muSun] = uplanet(cheap_departure, 3);
[kep2_min,muSun] = uplanet(cheap_intermediate, 6);

% transformation from radiants to degrees
kep1_min(3:6) = kep1_min(3:6)*180/pi;
kep2_min(3:6) = kep2_min(3:6)*180/pi;

% propagation of all the orbits
%orbit 1 earth
T1 = 2*pi * sqrt(kep1_min(1)^3/muSun);
tspan1 = linspace(0,T1,10000);
[r1,v1]=kepl_to_car(kep1_min(1),kep1_min(2),kep1_min(3),kep1_min(4),kep1_min(5),kep1_min(6),muSun);
y1 = [r1' v1'];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t1,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan1, y1, options);

%orbit 2 saturn
T2 = 2*pi * sqrt(kep2_min(1)^3/muSun);
tspan2 = linspace(0,T2,10000);
[r2,v2]=kepl_to_car(kep2_min(1),kep2_min(2),kep2_min(3),kep2_min(4),kep2_min(5),kep2_min(6),muSun);
y2 = [r2' v2'];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t2,Y2] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan2, y2, options);

%transfer orbit
ToF = (cheap_intermediate - cheap_departure)*24*3600;
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muSun, 0, 0, 0 );
TT = 2*pi * sqrt(A^3/muSun);
tspanT = linspace(0,TT,10000);
yT = [r1' VI];
options= odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[tT,YT1] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanT, yT, options);
tspanTOF = linspace(0,ToF,10000);
[tT,Yarc1] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanTOF, yT, options);
V_neg_s = VF;

%computation of deltav1_earth for the earth departure hyperbola
deltav1_earth = (VI'-v1);  %(V_satellite - V_earth) in the heliocentric frame, at cheap_departure

figure(3)
hold on
grid on
axis equal
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-');
plot3(Yarc1(:,1),Yarc1(:,2),Yarc1(:,3),'c-');
plot3(YT1(:,1),YT1(:,2),YT1(:,3),'c--');
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'bo');
legend('Earth orbit','Saturn orbit','transfer arc E-S','transfer orbit','Earth position','saturn position','Location','northeastoutside')
title('First heliocentric arc')

%% 5. heliocentric leg 2

% computation of algebraic ephemerides of the planets
[kep3_min,mass,M] = ephNEO(cheap_arrival,65);
[kep2_min,muSun] = uplanet(cheap_intermediate, 6);

% transformation from radiants to degrees
kep3_min(3:6) = kep3_min(3:6)*180/pi;
kep2_min(3:6) = kep2_min(3:6)*180/pi;

% propagation of all the orbits
%orbit 1
T3 = 2*pi * sqrt(kep3_min(1)^3/muSun);
tspan3 = linspace(0,T3,10000);
[r3,v3]=kepl_to_car(kep3_min(1),kep3_min(2),kep3_min(3),kep3_min(4),kep3_min(5),kep3_min(6),muSun);
y3 = [r3' v3'];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t3,Y3] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan3, y3, options);

%orbit 2 saturn
T2 = 2*pi * sqrt(kep2_min(1)^3/muSun);
tspan2 = linspace(0,T2,10000);
[r2,v2]=kepl_to_car(kep2_min(1),kep2_min(2),kep2_min(3),kep2_min(4),kep2_min(5),kep2_min(6),muSun);
y2 = [r2' v2'];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t2,Y2] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspan2, y2, options);

%transfer orbit
ToF2 = (cheap_arrival - cheap_intermediate)*24*3600;
[A,P,E,ERROR,VI2,VF2,TPAR,THETA] = lambertMR( r2, r3, ToF2, muSun, 0, 0, 0 );
TT = 2*pi * sqrt(A^3/muSun);
tspanT = linspace(0,TT,10000);
yT = [r2' VI2];
options= odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[tT,YT2] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanT, yT, options);
tspanTOF = linspace(0,ToF2,10000);
[tT,Yarc2] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanTOF, yT, options);
V_pos_s = VI2;
% heliocentric orbits plot

figure(4)
hold on
grid on
axis equal
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-');
plot3(Y3(:,1),Y3(:,2),Y3(:,3),'k-');
plot3(Yarc2(:,1),Yarc2(:,2),Yarc2(:,3),'g-');
plot3(YT2(:,1),YT2(:,2),YT2(:,3),'g--');
plot3(r2(1),r2(2),r2(3),'bo');
plot3(r3(1),r3(2),r3(3),'ko');
legend('Saturn orbit','Ganymed orbit','transfer arc S-G','transfer orbit S-G','Saturn position','Ganymed position','Location','northeastoutside')
title('Second heliocentric arc')

%% complete graph 

figure(5)
hold on
grid on
axis equal
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-');
plot3(Y3(:,1),Y3(:,2),Y3(:,3),'k-');
plot3(Yarc1(:,1),Yarc1(:,2),Yarc1(:,3),'c-');
plot3(YT1(:,1),YT1(:,2),YT1(:,3),'c--');
plot3(Yarc2(:,1),Yarc2(:,2),Yarc2(:,3),'g-');
plot3(YT2(:,1),YT2(:,2),YT2(:,3),'g--');
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'bo');
plot3(r3(1),r3(2),r3(3),'ko');

legend('Earth orbit','Saturn orbit','Ganymed orbit','transfer arc E-S','transfer orbit E-S',...
    'transfer arc S-G','transfer orbit S-G','Earth position',...
    'Ganymed position','saturn position','Location','northeastoutside')
title('Heliocentric trasfer representation')

%% 6. Earth sphere of influence

% circular orbit
% put a,e,...
R_earth = 6378.1;
a = 6778.1;
e = 0;
vcircular = sqrt(mu1/a);

% hyperbolic orbit
rp = a; %pericener of parabola is on the circular orbit
v_infinite_plus = norm(deltav1_earth); %comes from the first heliocentric leg

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

deltav_maneuver = abs(vp-vcircular);

% orbit propagation
%circular orbit
r_peri_c = rp.*[1;0;0]; 
r_peri_hyp = rp.*[cos(alpha);-sin(alpha);0];
yo = [r_peri_c' vcircular.*[0 1 0]]; %circular orbit state vector (pericenter)
y1 = [r_peri_hyp' vp.*[sin(alpha) cos(alpha) 0]]; %departing hyperbola state vector (pericenter)
%T = 10000000; 
tspano = linspace(10000,0,10000); %time span where the state vector is defined
tspan1 = linspace(0,1000,10000);
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,mu1), tspan1, y1, options);
[T,Yo] = ode113(@(t,y) twobodyproblem_ode(t,y,mu1), tspano, yo, options);
%Yo = Yo./R_earth;
%Y1 = Y1./R_earth;

figure(7)
grid on
hold on
axis equal
earth_sphere
plot3(Yo(:,1),Yo(:,2),Yo(:,3),'b-');
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
title('Earth SOI')

%% 7. Gravitational assist flyby

%V_pos_s from the Ganymed side (second helio leg, heliocentric frame), is
%defined above at line 304 via the lambert problem's solution

%V_neg_s from the Earth side (first helio leg, heliocentric frame), is
%defined above at line 251 via the lambert problem's solution

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

% DA AGGIUSTARE CON IL NUOVO FLYBY
%calculating approach and departure velocities
v_planet = v2; %velocity of saturn
v_inf_neg = V_neg_s'-v_planet;
v_inf_pos = V_pos_s'+v_planet;

%calculating the turning angle
v_neg = norm(v_inf_neg);
v_pos = norm(v_inf_pos);

turning_angle = acos(dot(v_inf_neg,v_inf_pos)/(v_pos*v_neg));

%calculating r_peri
turning_angle_function = @(r_peri) +(asin(1/(1+r_peri*(v_neg^2)/mu_planet)))+...
    (asin(1/(1+r_peri*(v_pos^2)/mu_planet)))-turning_angle/2

%definition of an initial guess
r_peri_guess = r_saturn*10; %Saturn's radius + random height, km
r_peri_fzero = fzero(turning_angle_function,r_peri_guess);
if abs(norm(r_peri_fzero)) < (r_saturn+100)
    fprintf('error, too low pericenter radius\n');
end
m_saturn = mu_planet/G;
m_sun = mu_sun/G;
R_saturn_mean = 9.5826*AU;
SOI = R_saturn_mean * (m_saturn/m_sun)^0.4;

if abs(r_peri_fzero) > (SOI)
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

alp = pi/2 - turning_angle/2;

%arcs plotting, planetocentric
r_peri = r_peri_fzero.*[sin(alp);-cos(alp);0]; 
yof = [r_peri' v_peri_in.*[cos(alp),sin(alp),0]]; %arrival state vector (pericenter)
y1f = [r_peri' v_peri_out.*[cos(alp),sin(alp),0]]; %departing state vector (pericenter)
%T = 10000000; 
tspano = linspace(1000000,0,100000); %time span where the state vector is defined
tspan1 = linspace(0,1000000,100000);
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y1f] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_planet), tspan1, y1f, options);
[T,Yof] = ode113(@(t,y) twobodyproblem_ode(t,y,mu_planet), tspano, yof, options);

figure(8)
grid on
hold on
axis equal
plot3(Yof(:,1),Yof(:,2),Yof(:,3),'b-');
plot3(Y1f(:,1),Y1f(:,2),Y1f(:,3),'r-');
plot3(0,0,0,'.k','MarkerSize',50);
legend('incoming','outcoming','saturn')

%%

r_peri_fzero=[0:1000:r_saturn*100]

% ADD TITLE AND LEGEND EVERYWHERE, LAST TWO SECTION TO ADJUST
a = -(2*asin(1./(1+r_peri_fzero.*(v_neg^2)./mu_planet)))/2-...
    (2*asin(1./(1+r_peri_fzero.*(v_pos^2)./mu_planet)))/2+turning_angle/2

figure(9)
plot(r_peri_fzero,a)