function pos_ecef = sat_kepler2ecef(a, e, i, omega, Omega, mu, t)
    
    % Constantes
    R_earth = 6378.137e3; % Radio ecuatorial medio

    % Posicion en ECI 
    pos_eci = sat_kepler2eci(a, e, i, omega, Omega, mu, t);

    %  matriz ECI to ECEF 
    theta = 7.2921159e-5 * t; % Omega tierra en rad/s
    R = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];

    % Transformacion from ECI to ECEF 
    pos_ecef = R * pos_eci;
    
end

% Comentarios: 

% Como la transformación es entre un sistema fijo ECI, y uno que rota a la
% velocidad de rotación de la tierra, es necesario tener en cuenta el
% tiempo que pasa entre muestra y muestra para que exista una correlación

