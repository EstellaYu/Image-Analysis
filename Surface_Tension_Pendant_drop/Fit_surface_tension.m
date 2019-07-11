function R = Fit_surface_tension(coeffs,Z,arcL)
%
% by W.D. Ristenpart
% March 31, 2006
%
% INPUTS
%  coeffs = 2 element vector, curvature and Bond number, which are to be optimized
%  Z = vector of smoothed Z coordinates
%  arcL = arclength of experimental curve (this sets the domain of integration)
%
% OUTPUT
%  R = numerically determined radii at each value of Z, by solving 'courbure' ODE
%  
%

% define parameters
Ro = coeffs(1);     % curvature  (pixels)
B = coeffs(2);      % Bond number (dimensionless)
L = arcL / Ro;      % length of integration (dimensionless)

% solve ODE for surface tension
[S,rzcoords] = ode45(@Courbure,[0 L],[0 0 0],[],B);

% rescale solution into pixels
znum = rzcoords(:,2) * Ro;
rnum = rzcoords(:,1) * Ro;

% interpolate solution to find a R value for each Z
R = interp1(znum,rnum,Z,'cubic');




% 
% 
% if length(a) == 1
%     Z = a(1)*R.^2;  
% elseif length(a) == 2
%     Z = a(2) + a(1)*R.^2;  
% elseif length(a) == 3
%     Z = a(2) + a(1)*(R-a(3)).^2;  
% end