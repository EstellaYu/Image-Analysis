function dr = Courbure(s,r,B)
%Here, we follow the method most ancient to the trade, which I can not for
%some reason, find from the force balance as set up in the book by deGennes
%et al. Anyway, the equation to solve is:
%
%   d\phi/ds = 2 - B x - sin(\phi)/y
%   dy/ds    = cos(\phi)
%   dx/ds    = sin(\phi)
%where x and y are the coordinates as you would expect, B is the bond
%number = (R_0)^2 \Delta \rho g / \gamma, where R_0 is the radius of curvature of the
%pendant drop at the apex (i.e. "1/b" in the previous program). \phi is the
%angle which the tangent to the drop at location (x,y) makes with the
%horizontal. ds is an infinitesimal length of the drop curve.

dr = zeros(3,1);% a column vector


if s==0
dr(1)=1;
dr(2)=1;
dr(3)=1;
else
dr(1)=cos(r(3));
dr(2)=sin(r(3));
dr(3)=2-sin(r(3))/r(1)-B*r(2);
end;