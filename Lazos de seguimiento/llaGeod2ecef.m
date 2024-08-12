function [r] = llaGeod2ecef(lla)
%LLAGEOD2ECEF Coordenadas ECEF de posición en latitud, longitud y altura
%   Convierte un conjunto de posiciones dadas en latitud, longitud y 
%	altitud geodésicas a sus coordenadas ECEF
% 
% ARGUMENTOS:
%	lla (3x1)	- Latitud geodésica, longitud y altitud [º,º,m]
% 
% DEVOLUCION:
%	p (KKx3)	- Coordenadas ECEF de las KK épocas expresadas en [m]

lat = lla(1);
lon = lla(2);
alt = lla(3);


a = 6378137;		% Semieje mayor del elipsoide WGS84
b = 6356752.3142;	% Semieje menor del elipsoide WGS84

e_2 = 1 - (b/a)^2;	% Excentricidad al cuadrado del elipsoide

coslon = cosd(lon);
coslat = cosd(lat);
sinlon = sind(lon);
sinlat = sind(lat);
tan2lat = (tand(lat)).^2;

x = (a.*coslon)./sqrt(1 + (1-e_2).*tan2lat) + alt.*coslon.*coslat;
y = (a.*sinlon)./sqrt(1 + (1-e_2).*tan2lat) + alt.*sinlon.*coslat;
z = (a.*(1-e_2).*sinlat)./sqrt(1 - e_2.*(sinlat.^2)) + alt.*sinlat;


r = [x; y; z];

end

