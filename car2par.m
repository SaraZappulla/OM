function [a,e,i,OM,om,th] = car2par(r, v,mu)

% Function used to compute keplerian parameters from the state vectors e
% and v expressed in Earth-centered equatorial frame

%INPUT
% r     position vector      [rx, ry, rz]'           [km]
% v     velocity vector      [vx, vy, vz]'           [km/s]

%OUTPUT
% a     semi major axes magnitude                           [km]
% e     eccentricity vector magnitude                       [adimensionale]
% i     inclination of orbital plane                        [Rad]
% OM    OMEGA = ascensione retta del nodo ascendente        [Rad]
% om    omega = anomalia del pericentro                     [Rad]
% th    theta = posizione angolare del corpo sull'orbita    [Rad]

%COSTANTI
Rt = 6378.1370; %[km]       Earth radius %astroConstants()
%mu = 398600.0;  %[km^3/s^2] Earth gravitarional parameter astroConstants()

%% 1.modulo dei vettori
R = norm(r);
V = norm(v); 

%% 2.calcolo a,e,i
a = 1/((2/R)-(V^2/mu));
h = cross(r,v);
H = norm(h);

%controllo che il modulo di h sia corretto
Hk = R*V;
if ((H-Hk)>10^3)
    fprintf('!error!\nmodulo di h non corrispondente %f %f', H, Hk);
end

e = (cross(v,h)/mu)-(r/R);
E = norm(e);
i = acos(h(3)/H);
n = cross([0 0 1]',h)/norm(cross([0 0 1]',h));

% define a vector that is:
eo = e; % if e not 0
% in case of circular orbit
if e == [0 0 0] 
    eo = n;
end

% h     momento angolare per unità di massa [km^2/s]
% H     modulo di h
% E     modulo di e
% n     asse dei nodi

%% 3.calcolo OM,om,th
if i == 0 % in case of equatorial orbit
    OM = 0;
else if n(2)>=0
    OM = acos(n(1));
else
    OM = 2*pi-acos(n(1));
end
end

%if e == [0 0 0]
 %   om = 0;
if (e(3)>=0)
    om = acos(dot(n,eo)/E);
else
    om = 2*pi-acos(dot(n,eo)/E);
end

Vr = dot(v,r)/R;
if Vr>0
    th = acos(dot(r,eo)/(R*E));
else
    th = 2*pi-acos(dot(r,eo)/(R*E));
end
end