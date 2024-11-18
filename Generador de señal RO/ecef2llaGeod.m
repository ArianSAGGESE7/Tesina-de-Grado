 function [lla] = ecef2llaGeod(r)
%ECEF2LLAGEOD Latitud, longitud y alturas geodéticas desde posición ECEF

% ARGUMENTOS:
%	r (3x1)		- Vector posición ECEF [m]
%	
% DEVOLUCION:
%	lla (3x1)	- Latitud, longitud y altitud [º], [º], [m]

x = r(1);
y = r(2);
z = r(3);

% Parámetros del elipsoide WGS84
a = 6378137;
b = 6356752.3142;
e = sqrt(1-(b/a)^2);
e2 = e^2;

lambda = atan2(y,x);

p = sqrt(x^2 + y^2);

phi_zm1 = 0;
phi = atan((z/p)/(1 - e2));

if (z ~= 0)
	while abs(phi_zm1 - phi) > 1e-9
		
		phi_zm1 = phi;
		
		y = sin(phi);
		N = a/sqrt(1 - e2*y^2);
		
		phi = atan((z/p)*(1 + e2*N*y/z));
	end
else
	y = 0;
end

if (z == 0) 
   y=0; 
end

h = p*cos(phi) + z*y - a*sqrt(1-e2*y^2);

lon = rad2deg(lambda);
lat = rad2deg(phi);
alt = h;

lla = [lat; lon; alt];
 
end

