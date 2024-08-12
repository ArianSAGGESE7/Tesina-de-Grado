function [x,v] = kepler2ecef(a, e, i, Omega, omega, nu)
    % Convertir parámetros de entrada de grados a radianes
    i = deg2rad(i);
    Omega = deg2rad(Omega);
    omega = deg2rad(omega);
    nu = deg2rad(nu);
    G = 6.67430e-11; % Constante de gravitación universal (m^3/kg/s^2)
    M = 5.972e24;    % Masa de la Tierra (kg)

    % Calcular la distancia radial
    % r = a * (1 - e^2) ./ (1 + e * cos(nu));

    % Coordenadas en el plano orbital
    % x_orbital = r .* cos(nu);
    % y_orbital = r .* sin(nu);

    % Rotar las coordenadas al sistema inercial
    % x_rotado = x_orbital .* (cos(omega) .* cos(Omega) - sin(omega) .* sin(Omega) .* cos(i)) - ...
    %            y_orbital .* (sin(omega) .* cos(Omega) + cos(omega) .* sin(Omega) .* cos(i));
    % y_rotado = x_orbital .* (cos(omega) .* sin(Omega) + sin(omega) .* cos(Omega) .* cos(i)) + ...
    %            y_orbital .* (cos(omega) .* cos(Omega) .* cos(i) - sin(omega) .* sin(Omega));
    % z_rotado = x_orbital .* (sin(omega) .* sin(i)) + y_orbital .* (cos(omega) .* sin(i));
    % 
    % % Devolver las coordenadas en metros
    % x = x_rotado ; % Convertir a metros
    % y = y_rotado ; % Convertir a metros
    % z = z_rotado ; % Convertir a metros
    r_orbital = a * (1 - e^2) ./ (1 + e * cos(nu));

    x_p = r_orbital.*cos(nu + omega);
    y_p = r_orbital.*sin(nu + omega);
    x_s = x_p.*cos(Omega)-y_p.*cos(i).*sin(Omega);
    y_s = x_p.*sin(Omega)+y_p.*cos(i).*cos(Omega);
    z_s = y_p.*sin(i);

    x=[x_s;y_s;z_s];


    gamma = atan(e*sin(nu)/(1+e*cos(nu)));

    v0_mod = sqrt(2*G*M/r_orbital-G*M/a);


    v_xs = v0_mod*(-1*sin(nu+omega-gamma)*cos(Omega) -1*cos(nu+omega-gamma)*cos(i)*sin(Omega));
    v_ys = v0_mod*(-1*sin(nu+omega-gamma)*sin(Omega) + cos(nu+omega-gamma)*cos(i)*cos(Omega));
    v_zs = v0_mod*(cos(nu+omega-gamma)*sin(i));
    
    v = [v_xs;v_ys;v_zs];

end
