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

%% Orbit data
clear
kep = [22594, 0.6847, 63.6246, 180,180,60];

%SRP data
p_sr = 4.5e-6; %at 1 AU, N/m^2
cr = 1.0; %reflection coefficient
Am_ratio = 5; %surface to mass ratio relative to the sun, m^2/kg 

%writing the functions
obl = astroConstants(8); %Earth's obliquity
mu = astroConstants(13); %earth
[r_spacecraft_earth,v_spacecraft_earth] = kepl_to_car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),mu);
r_spacecraft_sun = r_spacecraft_earth + [cos(obl);sin(obl);0]*astroConstants(2);
a_srp = p_sr*cr*Am_ratio*(astroConstants(2))^2/(norm(r_spacecraft_sun))^2;
acc_SRP = -a_srp.*r_spacecraft_sun./norm(r_spacecraft_sun)/1000;

%degrees to radians conversion, to be done beforehand
kep_rad = [kep(1) kep(2) 0 0 0 0];
kep_rad(3:6) = kep(3:6)*pi/180;
% kep in degrees, kep_rad in radiants

%rotation of acc_perturbed from RSW to carthesian ref. frame

%% 5. propagation of orbit with perturbation: expressed in cartesian coordinate and keplerian elements

% N of orbit considered
k = 100
% Gauss propagation
Tperiod = 2*pi * sqrt(kep(1)^3/mu);
tspan = linspace(0,Tperiod*k,1e5);
% s = state in this case written keplerian parameters
options = odeset('RelTol',1e-13,'AbsTol',1e-14); %specifications needed for ode113
[T,kep_propagate] = ode113(@(t,s) (EOM(t,s,mu,a_per_complete(t,s,acc_SRP))'), tspan, kep_rad, options);

% cartesian propagation only j2
[rr,vv]=kepl_to_car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),mu);
y1 = [rr;vv];
j2 = astroConstants(9); %earth gravitation perturbation
R = astroConstants(23); %radius of the Earth
[T,Y] = ode113(@(t,y) (twobodyproblem_perturbed_ode(t,y,mu,R,j2,acc_SRP)), tspan, y1, options);

% passing to keplerian parameters
aj2 = [];
ej2 = [];
ij2 = [];
OMj2 = [];
omj2 = [];
thetaj2 = [];
for j = 1:length(tspan)
    [aj2(j),ej2(j),ij2(j),OMj2(j),omj2(j),thetaj2(j)] = car_to_kepl(Y(j,1:3),Y(j,4:6),mu);
    %theta is always between 0 and 360 degrees
    if(thetaj2(j)>360)
        index = fix(thetaj2(j)/360);
        thetaj2(j) = thetaj2(j) - index*360;
    end
    if(OMj2(j)>360)
        index = fix(OMj2(j)/360);
        OMj2(j) = OMj2(j) - index*360;
    end
    if(omj2(j)>360)
        index = fix(omj2(j)/360);
        omj2(j) = omj2(j) - index*360;
    end
    
      if(abs(kep_propagate(j,6))>2*pi)
        index = fix(kep_propagate(j,6)/2*pi);
        kep_propagate(j,6) = kep_propagate(j,6)- index*2*pi;
      end
      if(kep_propagate(j,5)>2*pi)
        index = fix(kep_propagate(j,5)/2*pi);
        kep_propagate(j,5) = kep_propagate(j,5)- index*2*pi;
      end
      if(kep_propagate(j,4)>2*pi)
        index = fix(kep_propagate(j,4)/2*pi);
        kep_propagate(j,4) = kep_propagate(j,4)- index*2*pi;
      end

end

%% 6. Elements plotting

% a comparison plot
figure(1)
plot(tspan./Tperiod,kep_propagate(:,1));
hold on
plot(tspan./Tperiod,aj2);
legend('Gauss propagation','Carthesian propagation')
title('a variation in time')
xlabel('time[s]')
ylabel('a[km]')
grid on


%error plot
error = abs((kep_propagate(:,1)-aj2')./kep(1));
figure(2)
plot(tspan./Tperiod,error);
grid on
title('a methods error')
xlabel('time[s]')
ylabel('|a_G_a_u_s_s - a_C_a_r_t|')

%%
% e comparison plot
figure(3)
plot(tspan./Tperiod,kep_propagate(:,2));
hold on
plot(tspan./Tperiod,ej2);
legend('Gauss propagation','Carthesian propagation')
title('e variation in time')
xlabel('time[s]')
ylabel('e[-]')
grid on

%error plot
error = abs(kep_propagate(:,2)-ej2');
figure(4)
plot(tspan/Tperiod,error);
grid on
title('e methods error')
xlabel('time[s]')
ylabel('|e_G_a_u_s_s - e_C_a_r_t|[-]')


%% i comparison plot
figure(5)
plot(tspan/Tperiod,kep_propagate(:,3)*180/pi);
hold on
plot(tspan/Tperiod,ij2);
legend('Gauss propagation','Carthesian propagation')
title('i variation in time')
xlabel('time[s]')
ylabel('i[°]')
grid on
%error plot
error = abs(kep_propagate(:,3)*180/pi-ij2')/(360);
figure(6)
plot(tspan/Tperiod,error);
grid on
title('i methods error')
xlabel('time[s]')
ylabel('|i_G_a_u_s_s - i_C_a_r_t|')

%% OM comparison plot
figure(7)
plot(tspan/Tperiod,kep_propagate(:,4)*180/pi);
hold on
plot(tspan/Tperiod,OMj2);
legend('Gauss propagation','Carthesian propagation')
title('\Omega variation in time')
xlabel('time[s]')
ylabel('\Omega [°]')
grid on
%error plot
error = abs(kep_propagate(:,4)*180/pi-OMj2')/(360);
figure(8)
plot(tspan/Tperiod,error);
grid on
title('\Omega methods error')
xlabel('time[s]')
ylabel('|\Omega_G_a_u_s_s - \Omega_C_a_r_t|[-]')

%% om comparison plot
figure(9)
plot(tspan/Tperiod,kep_propagate(:,5)*180/pi);
hold on
plot(tspan/Tperiod,omj2);
legend('Gauss propagation','Carthesian propagation')
title('\omega variation in time')
xlabel('time[s]')
ylabel('\omega[°]')
grid on
%error plot
error = abs(kep_propagate(:,5)*180/pi-omj2')/(360);
figure(10)
plot(tspan/Tperiod,error);
grid on
title('\omega methods error')
xlabel('time[s]')
ylabel('|\omega_G_a_u_s_s - \omega_C_a_r_t|[-]')

%% theta comparison plot
figure(11)
plot(tspan/Tperiod,kep_propagate(:,6)*180/pi);
hold on
plot(tspan/Tperiod,thetaj2*180/pi);
legend('Gauss propagation','Carthesian propagation')
title('\theta variation in time')
xlabel('time[s]')
ylabel('\theta[°]')
grid on
%error plot
error = abs(kep_propagate(:,6)*180/pi-thetaj2')/(360);
figure(12)
plot(tspan/Tperiod,error);
grid on
title('\theta methods error')
xlabel('time[s]')
ylabel('|\theta_G_a_u_s_s - \theta_C_a_r_t|[-]')

%% 7. represent orbit evolution

% add a 3D plot of orbit
% just do a for-cycle to translate keplerian elements in rr
% kep_propagate are already made PA they are radiants 
% Do a plot3
kep_propagate_r = kep_propagate;
kep_propagate_d = kep_propagate;
kep_propagate_d(:,3:6) = kep_propagate(:,3:6).*180./pi;
for i = 1:length(kep_propagate)
    [yo(:,i) vo(:,i)] = kepl_to_car(kep_propagate_r(i,1),kep_propagate_r(i,2),kep_propagate_d(i,3),kep_propagate_d(i,4),kep_propagate_d(i,5),kep_propagate_d(i,6),mu);
end

figure(13)
plot3(Y(:,1),Y(:,2),Y(:,3),'r')
hold on
grid on
axis equal
plot3(yo(1,:),yo(2,:),yo(3,:),'b')
earth_sphere
legend('Carthesian propagation','Gauss perturbation')
title('Pertutbed orbits propagation')
xlabel('x[km]')
ylabel('y[km]')
zlabel('z[km]')



%% 8. filtering of high frequenies and plotting of both
N = 900; %or 900
points = 5000;
x = linspace(0,k*Tperiod,points);
[T,argument] = ode113(@(t,s) (EOM(t,s,mu,a_per_complete(t,s,acc_SRP))'), x, kep_rad, options);

% a graph
sec_a = movmean(argument(:,1),N);

figure(1)
plot(x./Tperiod,sec_a,'k');
legend('Gauss propagation','Carthesian propagation','Secular filter')
title('a variation in time')
xlabel('time[s]')
ylabel('a[km]')
grid on

%% e filtering
N = 900;
sec_e = movmean(argument(:,2),N);

figure(3)
plot(x./Tperiod,sec_e,'k');
legend('Gauss propagation','Carthesian propagation','Secular filter')
title('e variation in time')
xlabel('time[s]')
ylabel('e[-]')
grid on

%% i filtering
N = 100;
sec_i = movmean(argument(:,3),N);

figure(5)
plot(x./Tperiod,sec_i*180/pi,'k');
legend('Gauss propagation','Carthesian propagation','Secular filter')
title('i variation in time')
xlabel('time[s]')
ylabel('i[°]')
grid on

%% OM filtering
N = 100;
sec_OM = movmean(argument(:,4),N);

figure(7)
plot(x./Tperiod,sec_OM*180/pi,'k');
hold on; 

a_an = kep(1);
e_an = kep(2);
i_an = kep_rad(3);
OM_dot_analytic = -(3/2 * (j2.*sqrt(mu).*R.^2)/((1-e_an^2)^2*a_an^(7/2)))*cos(i_an);

[T,OM_dot_an] = ode113(@(t,s) (OM_dot_analytic), x, kep_rad(4), options);

figure(7)
plot(x./Tperiod,OM_dot_an*180/pi,'b');
legend('Gauss propagation','Carthesian propagation','Secular filter','Secular approximation')
title('\Omega variation in time')
xlabel('time[s]')
ylabel('\Omega [°]')
grid on

%% om filtering
N = 100;
sec_om = movmean(argument(:,5),N);

figure(9)
plot(x./Tperiod,sec_om*180/pi,'k');
hold on; 

om_dot_analytic = -(3/2 * (j2.*sqrt(mu).*R.^2)/((1-e_an^2)^2*a_an^(7/2)))*...
    (5/2*sin(i_an)^2-2);

[T,om_dot_an] = ode113(@(t,s) (om_dot_analytic), x, kep_rad(5), options);

figure(9)
plot(x./Tperiod,om_dot_an*180/pi,'b');
legend('Gauss propagation','Carthesian propagation','Secular filter','Secular approximation')
title('\omega variation in time')
xlabel('time[s]')
ylabel('\omega[°]')
grid on

%% theta filtering
N = 100;
% for j = 1:length(x)
%     if(argument(j,6)>2*pi)
%         index = fix(argument(j,6)/2*pi);
%         argument(j,6) = argument(j,6)- index*2*pi;
%     end
% end
sec_theta = movmean(argument(:,6),N);

figure(11)
plot(x./Tperiod,sec_theta*180/pi,'k');
legend('Gauss propagation','Carthesian propagation','Secular filter')
title('\theta variation in time')
xlabel('time[s]')
ylabel('\theta[°]')
grid on

%% 9. select an object and download the orbital elements - propagate its
% orbit using our model - compare real and our model