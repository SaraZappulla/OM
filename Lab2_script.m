%Second laboratory (12 Oct. 2022)

%%
%Exercise 1 - Representing the groundtracks

close all;
clear;
clc;

%Initial conditions and constant terrestrial parameters
%rr = [-4578.219;-801.084;-7929.708]; %km
%vv = [0.8; -6.037; 1.385]; %km/s
a = 7171.01; %km
e = 0;
in = 98;
OM = 0;
om = 40;
theta = 0;

%a = 26600;
%e = 0.74;
%in = 63.4;
%OM = 50;
%om = 280;
%theta = 0;
thetag_zero = 0; %degrees
time_vector = (linspace(0,80000,9000))'; 
omega_e = 15.04/3600; %deg/s, depends on time_vector spacing
mu_earth = astroConstants(13);

%calculation of the orbital period
%[a,e,in,OM,om,theta] = car_to_kepl(rr,vv,mu_earth);
T = 2*pi*sqrt((a^3)/mu_earth);
time_vector = (linspace(0,3*T,100000))'; 

%groundtrack function
[alpha,delta,lon,lat] = groundtrack(mu_earth,thetag_zero,time_vector,omega_e,a,e,in,OM,om,theta);
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
plot(lon180,lat180,'LineStyle','none','Marker','.','Color','g','MarkerSize',2)

%%
%Exercise 2 - Repeating ground tracks

k=15;
m=1;
[semimaj_axis] = repeating_groundtrack_semimajor_axis(mu_earth,omega_e,m,k);
T = 2*pi*sqrt((a^3)/mu_earth);
time_vector = (linspace(0,3*T,100000))'; 
[alpha,delta,lon,lat] = groundtrack(mu_earth,thetag_zero,time_vector,omega_e,semimaj_axis,e,in,OM,om,theta);
lon180 = wrapTo180(lon);
lat180 = wrapTo180(lat);

%singularity plotting conditions on latitude and longitude
lat(end+1) = inf;
lon(end+1) = inf;

figure(2)
xlin = [-180 180];
ylin = [-90 90];
I = imread("Earth for groundtracks.jpg","jpg");
h = image('Xdata',xlin,'Ydata',ylin,'Cdata',I);
hold on
plot(lon180,lat180,'LineStyle','none','Marker','.','Color','r','MarkerSize',2)

%%
%Exercise 3 - Ground track from ephemerides table

%Loading the ISS's ephemerides matrix file from the folder
%NOTE: change the data to match slides' data (5/10/22 - 6/10/22)
load("Ephemerides_for_ISS.mat");

%Imported file description: 
%eccentricity in the first column
%inclination in the second column
%OM in the third column
%om in the fourth column
%theta in the fifth column
%a in the sixth column
%
%Every row "m" is a point of the orbit sampled at a time t = m*n_row (minutes)

%Definition of empty rr and vv vectors for each point (1441x3 size each)
z = EphemeridesforISS(:,1);
r = zeros(size(z,1),3);
v = zeros(size(z,1),3);


time_vector = (linspace(0,size(z,1),size(z,1)))';

for k=1:size(z,1)
    %parameters transport for shorter conversion formula
    e = EphemeridesforISS(k,1);
    in = EphemeridesforISS(k,2);
    OM = EphemeridesforISS(k,3);
    om = EphemeridesforISS(k,4);
    theta = EphemeridesforISS(k,5);
    a = EphemeridesforISS(k,6);

    %Conversion to cartesian coordinates
    [r(k,:),v(k,:)] = kepl_to_car(a,e,in,OM,om,theta,mu_earth);

    %Calculating latitude and longitude
    [~,~,lon,lat] = groundtrack(mu_earth,thetag_zero,time_vector,omega_e*60,a,e,in,OM,om,theta);
    lon_out(k) = lon(k);
    lat_out(k) = lat(k);
end

lon180 = wrapTo180(lon_out);
lat180 = wrapTo180(lat_out);

%add kelpersolver of ISS
figure(3)
xlin = [-180 180];
ylin = [-90 90];
I = imread("Earth for groundtracks.jpg","jpg");
h = image(xlin,ylin,I);
hold on
plot(lon180,lat180,'LineStyle','none','Marker','.','Color','r','MarkerSize',15)
plot(lon180,lat180,'LineStyle','none','Marker','.','Color','r','MarkerSize',15)


