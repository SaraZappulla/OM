function [E] = KeplerEqSolver(tspan,e,a,mu)

%
%
%
n = sqrt(mu/((a)^3));
%if nargin<7 & nargin > 5
 %   if Eo <=0 || Eo > 1e-4
toll = 1e-6;
   % else 
    %    toll = Eo;
Eo = n*tspan + e*sin(n*tspan)/(1-sin(n*tspan+e)+sin(n*tspan));   
   % end
%end
E = [];
for i = 1:length(tspan)
    E(i) = fzero(@(E) n*tspan(i)-E+e*sin(E),Eo);
end