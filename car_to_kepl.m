function[a,e,i,OM,om,theta] = car_to_kepl(rr,vv,mu)

% FUNCTION DESCRIPTION
% this function provide orbit parameters in perifocal frame given in input
% carthesia variable

% INPUT


% OUTPUT


%AUTHORS
%angle outputs are in degrees!
r=norm(rr);
v=norm(vv);
E=0.5*v^2-mu/r;
a=-mu/(2*E);
hh=cross(rr,vv);
h=norm(hh);
i=acos(hh(3)/h)*(180/pi);
NN=cross([0;0;1],hh);
N=norm(NN);
ee=(1/mu).*((v^2-(mu/r)).*rr-dot(rr,vv).*vv);
e=norm(ee);
if (dot(vv,rr)>=0)
    theta=acos((dot(rr,ee))/(r*e))*(180/pi);
else
    theta=360-acos((dot(rr,ee))/(r*e))*(180/pi);
end 
if (dot(ee,[0;0;1])>=0)
    om=acos((dot(NN,ee))/(N*e))*(180/pi);
else
    om=360-acos((dot(NN,ee))/(N*e))*(180/pi);
end
if (dot(NN,[0;1;0])>=0)
    OM=acos(NN(1,1)/N)*(180/pi);
else
    OM=360-acos(NN(1,1)/N)*(180/pi);
end
end