%% Exercise 1 part 1 - Gaussian CORRECT, DO NOT CHANGE
clear
%SRP data
p_sr = 4.5e-6; %at 1 AU, N/m^2
cr = 1.0; %reflection coefficient
Am_ratio = 5; %surface to mass ratio relative to the sun, m^2/kg 

%writing the functions
obl = astroConstants(8); %Earth's obliquity
mu = astroConstants(13); %earth
kep = [7571,0.01,87.9,180,180,0]; %angles are already in degrees
[r_spacecraft_earth,v_spacecraft_earth] = kepl_to_car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),mu);
r_spacecraft_sun = r_spacecraft_earth + [cos(obl);sin(obl);0]*astroConstants(2);
a_srp = p_sr*cr*Am_ratio*(astroConstants(2))^2/(norm(r_spacecraft_sun))^2;
acc_SRP = -a_srp.*r_spacecraft_sun./norm(r_spacecraft_sun)./1000;

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

%rotation of acc_perturbed from RSW to carthesian ref. frame

%% Exercise 1 part 2 - Carthesian CORRECT

% comparison with two body problem perturbed
[rr,vv]=kepl_to_car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),mu);
y1 = [rr;vv];
j2 = astroConstants(9);
mu = astroConstants(13);
R = astroConstants(23);
[T,Y] = ode113(@(t,y) (twobodyproblem_perturbed_ode(t,y,mu,R,j2,acc_SRP)), tspan, y1, options);

a = [];
e = [];
i = [];
OM = [];
om = [];
theta = [];
for j = 1:length(tspan)
    [a(j),e(j),i(j),OM(j),om(j),theta(j)] = car_to_kepl(Y(j,1:3),Y(j,4:6),mu);
    %theta is always between 0 and 360 degrees
    if(theta(j)>360)
        index = fix(theta(j)/360);
        theta(j) = theta(j) - index*360;
    end
end


%% e comparison plot
figure(3)
plot(tspan./Tperiod,kep_propagate(:,2));
hold on
plot(tspan./Tperiod,e);
legend('Gauss','two body perturbed')
grid on

%error plot
error = abs(kep_propagate(:,2)-e');
figure(4)
plot(tspan/Tperiod,error);


%% a comparison plot
figure(1)
plot(tspan/Tperiod,kep_propagate(:,1));
hold on
plot(tspan/Tperiod,a);
legend('Gauss','carthesian')

%error plot
error = abs((kep_propagate(:,1)-a')./kep(1));
figure(2)
plot(tspan/Tperiod,error);
%% i comparison plot
figure(5)
plot(tspan/Tperiod,kep_propagate(:,3)*180/pi);
hold on
plot(tspan/Tperiod,i);
legend('Gauss','carthesian')

%error plot
error = abs(kep_propagate(:,3)*180/pi-i')/(360);
figure(6)
plot(tspan/Tperiod,error);

%% OM comparison plot
figure(7)
plot(tspan/Tperiod,kep_propagate(:,4)*180/pi);
hold on
plot(tspan/Tperiod,OM);
legend('Gauss','two body perturbed')

%error plot
error = abs(kep_propagate(:,4)*180/pi-OM')/(360);
figure(8)
plot(tspan/Tperiod,error);

%% om comparison plot
figure(9)
plot(tspan/Tperiod,kep_propagate(:,5)*180/pi);
hold on
plot(tspan/Tperiod,om);
legend('Gauss','two body perturbed')

%error plot
error = abs(kep_propagate(:,5)*180/pi-om')/(360);
figure(10)
plot(tspan/Tperiod,error);

%% theta comparison plot
figure(11)
plot(tspan/Tperiod,kep_propagate(:,6)*180/pi);
hold on
plot(tspan/Tperiod,theta*180/pi);
legend('Gauss','two body perturbed')

%error plot
error = abs(kep_propagate(:,6)*180/pi-theta')/(360);
figure(12)
plot(tspan/Tperiod,error);

%%
% fun = @(x) x + 5*cos(x) + 2*sin(4*x);
%x = linspace(0,8*pi,360*4+1)
%y = fun(x)

%N1 = 90
%ymean1 = movmean(y,N1)

%N2 = 360
%ymean2 = movmean(y,N2)
%% Exercise 1 part 3: a filtering

N = 200;
x = linspace(0,100*Tperiod,N);
%[T,argument] = ode113(@(t,s) (EOM(t,s,mu,a_per(t,s))'), x, kep_rad, options);

%sec_a = movmean(argument(:,1),N);

%200 filter the major oscillation
%2 doen't work
%4 still little oscillation
%3 is the maximum
N = 200;
z = linspace(0,100*Tperiod,5000);
[T,argument1] = ode113(@(t,s) (EOM(t,s,mu,a_per_complete(t,s,acc_SRP))'), z, kep_rad, options);
sec_b = movmean(argument1(:,1),N);

figure(1)
plot(z./Tperiod,sec_b);

%% e filtering

N = 100;
punti = 5000;
x = linspace(0,100*Tperiod,punti);
[T,argument] = ode113(@(t,s) (EOM(t,s,mu,a_per(t,s))'), x, kep_rad, options);

% sec_a = movmean(argument(:,2),N);

%300 filters the major oscillation
%100 filters once a period
%2 filters twice in 100 but we still have oscillations
%8 filters this last oscillation
%4 it also works without oscillations
% 
% N = 4;
% x = linspace(0,100*Tperiod,N);
% [T,argument] = ode113(@(t,s) (EOM(t,s,mu,a_per(t,s))'), x, kep_rad, options);

sec_b = movmean(argument(:,2),N);

figure(3)
plot(x./Tperiod,sec_b);

%% i filtering

% N = 300;
% x = linspace(0,100*Tperiod,N);
% [T,argument] = ode113(@(t,s) (EOM(t,s,mu,a_per(t,s))'), x, kep_rad, options);
% 
% sec_a = movmean(argument(:,3),N);

%300 filters the major oscillation
%100 filters once a period
%2 filters twice in 100 but we still have oscillations
%8 filters this last oscillation
%4 it also works without oscillations

N = 3;
x = linspace(0,100*Tperiod,N);
[T,argument] = ode113(@(t,s) (EOM(t,s,mu,a_per(t,s))'), x, kep_rad, options);

sec_b = movmean(argument(:,3),N);

figure(5)
plot(x/Tperiod,sec_b*180/pi,'k');

%% OM filtering
N = 500;
x = linspace(0,100*Tperiod,N);
[T,argument] = ode113(@(t,s) (EOM(t,s,mu,a_per(t,s))'), x, kep_rad, options);
sec_b = movmean(argument(:,4),N);

figure(7)
plot(x/Tperiod,sec_b*180/pi,'k');
hold on; 

a_an = kep(1);
e_an = kep(2);
i_an = kep_rad(3);
OM_dot_analytic = -(3/2 * (j2.*sqrt(mu).*R.^2)/((1-e_an^2)^2*a_an^(7/2)))*cos(i_an);

[T,OM_dot_an] = ode113(@(t,s) (OM_dot_analytic), x, kep_rad(4), options);

figure(7)
plot(x/Tperiod,OM_dot_an*180/pi,'b');

%% om filtering
N = 3;
x = linspace(0,100*Tperiod,N);
[T,argument] = ode113(@(t,s) (EOM(t,s,mu,a_per(t,s))'), x, kep_rad, options);
sec_b = movmean(argument(:,5),N);

figure(9)
plot(x/Tperiod,sec_b*180/pi,'k');
hold on; 

om_dot_analytic = -(3/2 * (j2.*sqrt(mu).*R.^2)/((1-e_an^2)^2*a_an^(7/2)))*...
    (5/2*sin(i_an)^2-2);

[T,om_dot_an] = ode113(@(t,s) (om_dot_analytic), x, kep_rad(5), options);

figure(9)
plot(x/Tperiod,om_dot_an*180/pi,'b');

%% 
