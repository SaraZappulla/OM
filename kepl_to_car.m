function[rr,vv]=kepl_to_car(a,e,i,OM,om,theta,mu)
%convertitore da parametri kepleriani a coordinate cartesiane. Gli angoli
%devono essere inseriti espressi in gradi, gli altri parametri orbitali
%sono scalari.
i=i.*(pi/180);
om=om.*(pi/180);
OM=OM.*(pi/180);
theta=theta.*(pi/180);

p=a.*(1-e.^2);
r=p./(1+e.*cos(theta));

vr=sqrt(mu./p).*e.*sin(theta);
va=sqrt(mu./p)*(1+e.*cos(theta));
rrpf=[r.*cos(theta);r.*sin(theta);0];
vvpf=[vr.*cos(theta)-va.*sin(theta);vr.*sin(theta)+va.*cos(theta);0];

MR=[cos(om) sin(om) 0;-sin(om) cos(om) 0;0 0 1]*[1 0 0;0 cos(i) sin(i);...
    0 -sin(i) cos(i)]*[cos(OM) sin(OM) 0;-sin(OM) cos(OM) 0;0 0 1];
MR=MR';
rr=MR*rrpf;
vv=MR*vvpf;
end