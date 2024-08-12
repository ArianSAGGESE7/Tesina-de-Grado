%% Pruebas del propagador


close all;clear all;clc

% algunas constantes de utilidad
G = 6.67430e-11; % Constante de gravitación universal (m^3/kg/s^2)
M = 5.972e24;    % Masa de la Tierra (kg)
R = 6371000;     % Radio de la Tierra (m)


% Datos Keplerianos de la orbita

%GPS

% a = 20000e3; % Semieje mayor en kilómetros
% e = 0.0; % Excentricidad
% i = 45; % Inclinación en radianes
% Omega = -60; % Longitud del nodo ascendente en radianes
% omega = 0; % Argumento del perigeo
% nu = 150;

%LEO

a = R+600e3; % Semieje mayor en kilómetros
e = 0.0; % Excentricidad
i = 97.77; % Inclinación en radianes
Omega = 0; % Longitud del nodo ascendente en radianes
omega = 0; % Argumento del perigeo
nu = 0;



% Otros parámetros

Periodo=2*pi*sqrt((a)^3/(M*G));
resolucion = 1 ; % en segundos
vueltas = 2    ;

%% KEPLER -> ECI -> x0, v0 -> resolución de EDO
% Esta función trabaja en coordenadas ECI, se le da al propagador una
% posición y una velocidad inicial y a partir de ello resuelve la EDO

[x,v] = Propagador_de_orbita(a, e, i, Omega, Omega, nu,Periodo,resolucion,vueltas);

%% KEPLER -> ECI -> ECEF -> (Para una dada especificación del tiempo)

% A diferencia del la anterior función esta resuelve la posición en ECEF de
% cada punto temporal que se le dá como referencia (recordar que tiene toe = tau = perigeo = 0)
% el contador va hasta que se cumple un período de orbita
% La resolución es de un segundo es decir que se resuelve la posición del
% satélite cada un segundo, esto puede variar cambiando "resolución" a su
% vez podemos variar la cantidad de vueltas que damos

inc = 45; % Inclinación en radianes

for i=1:(floor(Periodo))*vueltas+resolucion
    pos_ecef = sat_kepler2ecef(a, e, deg2rad(inc), deg2rad(omega), deg2rad(Omega), G*M, (i-1)*resolucion);
    pos_pos(i,:) = pos_ecef';
end
tierra_color()
comet3(pos_pos(:, 1), pos_pos(:, 2), pos_pos(:, 3))

%% ECEF -> lat,lon,h -> Ground track

% Cálculo del Ground track sobre la tierra a partir de encontrar la latitud
% longitud y altura

for i=1:length(pos_pos)
    [lla1] = ecef2llaGeod(pos_pos(i,:)); 
    lla(i,:) = lla1;
end
figure;

% geoplot(lla(:,1)',lla(:,2)','-*')
worldmap('World')
load coastlines
plotm(coastlat,coastlon)
hold on
plotm(lla(:,1)',lla(:,2)'+90,'-o','Color','k')

%%
figure 
geoplot(lla(:,1)',lla(:,2)'+90,'-*')