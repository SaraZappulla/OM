% Take Lambert function
clc 
clear
%% exercise 1 computing lambert function
muEarth = astroConstants(13);      % Sun's gravitational parameter [km^3/s^2];
ToF = 15*3600+6*60+40;         % Time in [s];
r1 = [-21800,37900,0];       % Initial position vector [km]
r2 = [27300,27700,0];     % Final position vector [km]
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muEarth, 0, 0, 0 );

% EX1 part 2 computing the all orbit
yo = [r1 VI];
y1 = [r2 VF];
T = 2*pi * sqrt(A^3/muEarth);
tspan = linspace(0,T,10000); %%time span where the state vector is defined
%[dy] = twobodyproblem_ode(tspan, y1,muEarth)
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,Y] = ode113(@(t,y) twobodyproblem_ode(t,y,muEarth), tspan, y1, options);

figure(1)
earth_sphere;
axis equal
grid on
hold on;
plot3(Y(:,1),Y(:,2),Y(:,3),'k-');
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'bo');

%% exercise 2
muEarth = astroConstants(13);
kep1 = [12500; 0; 0; 0; 0; 120];
kep2 = [9500; 0.3; 0; 0; 0; 250];
ToF = 3300; 
[r1,v1]=kepl_to_car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),muEarth);
[r2,v2]=kepl_to_car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muEarth);
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muEarth, 0, 0, 0 );
kep_t = [A,E,0,0,0,THETA]; %transfer orbit parameters
deltav1 = VI'-v1;
deltav2 = v2-VF';
deltavtot = norm(deltav1)+norm(deltav2);
ndeltav1 = norm(deltav1);
ndeltav2 = norm(deltav2);

%computing all the orbits
%orbit 1
T1 = 2*pi * sqrt(kep1(1)^3/muEarth);
tspan1 = linspace(0,T1,10000);
y1 = [r1 v1];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t1,Y1] = ode113(@(t,y) twobodyproblem_ode(t,y,muEarth), tspan1, y1, options);

%orbit 2
T2 = 2*pi * sqrt(kep2(1)^3/muEarth);
tspan2 = linspace(0,T2,10000);
y2 = [r2 v2];
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[t2,Y2] = ode113(@(t,y) twobodyproblem_ode(t,y,muEarth), tspan2, y2, options);

%transfer orbit
TT = 2*pi * sqrt(A^3/muEarth);
tspanT = linspace(0,TT,10000);
yT = [r1' VI];
options= odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[tT,YT] = ode113(@(t,y) twobodyproblem_ode(t,y,muEarth), tspanT, yT, options);
tspanTOF = linspace(0,ToF,10000);
[tT,Yarc] = ode113(@(t,y) twobodyproblem_ode(t,y,muEarth), tspanTOF, yT, options);


figure(2)
earth_sphere;
axis equal
grid on
hold on;
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-');
plot3(Yarc(:,1),Yarc(:,2),Yarc(:,3),'c-');
plot3(YT(:,1),YT(:,2),YT(:,3),'c--');
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'bo');
legend('initial orbit','final orbit','transfer arc','transfer orbit','initial transfer point','final transfer point','Location','northeastoutside')

%% exercise 3
dep_time1 = date2mjd2000([2003,4,1,12,0,0]);
dep_time2 = date2mjd2000([2003,8,1,12,0,0]);
arr_time1 = date2mjd2000([2003,9,1,12,0,0]);
arr_time2 = date2mjd2000([2004,3,1,12,0,0]);
time_length1_days = (dep_time2 - dep_time1); % length of vector in days
time_length1_sec = time_length1_days * 24 * 3600; % length of vector in seconds
time_length2_days = (arr_time2 - arr_time1); % vectors of days
time_length2_sec = time_length2_days * 24 *3600; % length of vector in seconds
departing_window = linspace(dep_time1,dep_time2,time_length1_days); %dep window expressed in seconds
arrival_window = linspace(arr_time1,arr_time2,time_length2_days); %arr window expressed in seconds
mu1 = astroConstants(13);
mu2 = astroConstants(14);


% for done for each second of the starting window and arrival window
time_line = [];
deltavtot = [];
for i = 1:time_length1_days
    for j = 1:time_length2_days
        t1 = departing_window(i)*3600*24; %seconds
        t2 = arrival_window(j)*3600*24; %seconds
        [kep1,ksun] = uplanet(departing_window(i), 3);
        [kep2,ksun] = uplanet(arrival_window(j), 4);
        kep1(6) = kep1(6)*180/pi;
        kep2(6) = kep2(6)*180/pi;
        kep1(3) = kep1(3)*180/pi;
        kep2(3) = kep2(3)*180/pi;
        kep1(4) = kep1(4)*180/pi;
        kep2(4) = kep2(4)*180/pi;
        kep1(5) = kep1(5)*180/pi;
        kep2(5) = kep2(5)*180/pi;
        deltavtot(i,j) = deltavtot_computation(t1,t2,kep1,kep2,mu1,mu2,ksun);
        time_line(i,j) = (arrival_window(j) - departing_window(i));
    end
end

%
time_difference_6 = [];
time_difference_12 = [];
time_difference_18 = [];
time_difference_24 = [];
time_difference_30 = [];
k_6 = 1;
k_12 = 1;
k_18 = 1;
k_24 = 1;
k_30 = 1;
for i = 1:time_length1_days
    for j = 1:time_length2_days
        if ((time_line(i,j) > 59.9) && (time_line(i,j) < 61))
            time_difference_6(1,k_6) = departing_window(i);
            time_difference_6(2,k_6) = arrival_window(j);
            k_6 = k_6+1;
        end
        if ((time_line(i,j) > 119.9) && (time_line(i,j) < 121))
            time_difference_12(1,k_12) = departing_window(i);
            time_difference_12(2,k_12) = arrival_window(j);
            k_12 = k_12+1;
        end
         if ((time_line(i,j) > 179.999) && (time_line(i,j) < 181))
           time_difference_18(1,k_18) = departing_window(i);
            time_difference_18(2,k_18) = arrival_window(j);
            k_18 = k_18+1;
         end
         if ((time_line(i,j) > 239.9) && (time_line(i,j) < 241))
            time_difference_24(1,k_24) = departing_window(i);
            time_difference_24(2,k_24) = arrival_window(j);
            k_24 = k_24+1;
         end
         if ((time_line(i,j) > 299.9) && (time_line(i,j) < 301))
            time_difference_30(1,k_30) = departing_window(i);
            time_difference_30(2,k_30) = arrival_window(j);
            k_30 = k_30+1;
        end
    end
end

%% optimisation graph

% to obtain dates inside the graph
dep_date = [];
arr_date = [];
for i=1:time_length1_days
    dep_date(i,:) = mjd20002date(departing_window(i));
end
for i=1:time_length2_days
    arr_date(i,:) = mjd20002date(arrival_window(i));
end

A = num2str(dep_date(:,1:3)); %yymmdd
B = num2str(arr_date(:,1:3));

figure(3)
contour(departing_window,arrival_window,deltavtot',linspace(5.5,14,25));
c = colorbar;
hold on
grid on
plot(cheap_departure,cheap_arrival,'or') %cheap transfer
plot(time_difference_6(1,:),time_difference_6(2,:),'-k');
plot(time_difference_12(1,:),time_difference_12(2,:),'-k');
plot(time_difference_18(1,:),time_difference_18(2,:),'-k');
plot(time_difference_24(1,:),time_difference_24(2,:),'-k');
plot(time_difference_30(1,:),time_difference_30(2,:),'-k');
c.Label.String = '\DeltaV[km/s]';
xticklabels(A);
yticklabels(B);
%contour(departing_window,arrival_window,time_line',[60,120,180,240,300],'k');
%colorbar(z,'off')
xlabel('Departing Date')
ylabel('Arrival Date')

%% minimum delta v transfer

mindelta = min(min(deltavtot))

for i=1:time_length1_days
    for j = 1:time_length2_days
        if deltavtot(i,j) == mindelta
            k=i
            c=j
        end
    end
end

cheap_departure = departing_window(k);
cheap_arrival = arrival_window(c);
cheap_departure_date = mjd20002date(cheap_departure);
cheap_arrival_date = mjd20002date(cheap_arrival);

% computation of algebraic ephemerides of the planets
[kep1_min,muSun] = uplanet(cheap_departure, 3);
[kep2_min,muSun] = uplanet(cheap_arrival, 4);

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
ToF = (cheap_arrival - cheap_departure)*24*3600;
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muSun, 0, 0, 0 );
TT = 2*pi * sqrt(A^3/muSun);
tspanT = linspace(0,TT,10000);
yT = [r1' VI];
options= odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[tT,YT] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanT, yT, options);
tspanTOF = linspace(0,ToF,10000);
[tT,Yarc] = ode113(@(t,y) twobodyproblem_ode(t,y,muSun), tspanTOF, yT, options);

%% orbits plot

figure(4)
hold on
grid on
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'r-');
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b-');
plot3(Yarc(:,1),Yarc(:,2),Yarc(:,3),'c-');
plot3(YT(:,1),YT(:,2),YT(:,3),'c--');
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'bo');
legend('initial orbit','final orbit','transfer arc','transfer orbit','initial transfer point','final transfer point','Location','northeastoutside')
