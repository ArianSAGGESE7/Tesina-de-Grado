%% Script para cálculo de ángulos 

% La idea de este script es tratar de encontrar los ángulos de interés para
% el patrón de radiación de la antena por medio de una transformación de
% espacios/ o de ejes coordenados.

clc;clear;close all;

% Llevamos todo al rededor del 0 

% El punto final le resto la posición de su correspondiente LEO (nuestro USAT)

rGPS = load('POSICIONES_PRN.mat');
rLEO = load('POSICIONES_USAT.mat');
velLEO = load('velocidad_usat.mat');
Posicion_final_ray = load('ultima_posicion.mat');
Exitosos = load('exitosos.mat');

rGPS=rGPS.EVENTO_PRN_fino;
Posicion_LEO = rLEO.EVENTO_USAT_fino;
Velocidad_LEO = velLEO.EVENTO_PRN_vel;
Posicion_final_ray = Posicion_final_ray.ultima_posicion;
Exitosos=Exitosos.exitosos; % contiene cuales de las posiciones de USAT 
% son las que se pudieron resolver con los scripts de RO
contador =1;
for i = Exitosos

    % Siempre vamos a estar centrados en 0 pero al punto final el vamos a
    % restar un distinto punto de ubicación del USAT
    Punto = Posicion_final_ray(:,contador) - Posicion_LEO(:,i); % llevo al cero conservando el ángulo 
    [T_org2usat] = transformacion_USAT(Posicion_LEO(:,i),Velocidad_LEO(:,i));
    Punto_transformado = T_org2usat*Punto;
 
    % Ahora nuestro eje x es la dirección que apunta a la tierra
    % nuestro eje y es la dirección instantanea de velocidad
    % y nuestro eje z es la terna complementaria
    
    % me interesa ver el ángulo sin componente y (elevación)
    % me interesa ver el ángulo sin componente z (azimut)

    % Calculamos cada uno de esos ángulos normalizando ese vector,
    % conservando dichos ángulos
    Punto_transformado = Punto_transformado/norm(Punto_transformado);
    ang_el = atan(Punto_transformado(3)/Punto_transformado(1));
    ang_azimut = atan(Punto_transformado(2)/Punto_transformado(1));

end