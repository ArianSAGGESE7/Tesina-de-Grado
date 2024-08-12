function pos_eci = sat_kepler2eci(a, e, i, omega, Omega, mu, t)

    % Movimiento medio (freq)
    n = sqrt(mu / a^3);

    % Anomalía media
    M = n * (t - 0); % Assuming epoch at t = 0

    % Anomalia Excentrica (resuelta iterativamente)
    E = solve_keplers_equation(M, e);

    % Anomalia verdadera
    nu = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

    % Raio orbital
    r = a * (1 - e * cos(E));

    % Posicion en el plano orbital (perifocal)
    x_orbital = r * cos(nu);
    y_orbital = r * sin(nu);
    pos_orbital = [x_orbital; y_orbital; 0];

    % Transformación a coordenadas ECI
    Q = [cos(omega), -sin(omega) * cos(i), sin(omega) * sin(i);
         sin(omega), cos(omega) * cos(i), -cos(omega) * sin(i);
         0, sin(i), cos(i)];

    % Rotate position to ECI 
    pos_eci = Q * pos_orbital;

    % Rotación debido a la ascension del nodo Omega
    R = [cos(Omega), -sin(Omega), 0;
         sin(Omega), cos(Omega), 0;
         0, 0, 1];
    pos_eci = R * pos_eci;
end
