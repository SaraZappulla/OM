function [deltavtot] = deltavtot_computation(t1,t2,kep1,kep2,mu1,mu2,muSun)

%INPUT
%kep1 = [a,e,i,OM,om,theta1] with theta1 dummy value in the script but 
% defined in the function

%OUTPUT
%
%

ToF = t2-t1; 

[r1,v1]=kepl_to_car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),muSun);
[r2,v2]=kepl_to_car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muSun);
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR( r1, r2, ToF, muSun, 0,0,0, 0);

deltav1 = VI'-v1;
deltav2 = v2-VF';
deltavtot = abs(norm(deltav1)+norm(deltav2));

vcircle = sqrt(mu1/kep1(1));
vinf=sqrt(2)*vcircle;

if deltav1<=vinf
    deltavtot = 0;
end

return