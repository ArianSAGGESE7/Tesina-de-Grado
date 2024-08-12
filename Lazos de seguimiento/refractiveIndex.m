% Funcion para calcular el indice de refracci�n n seg�n el punto en el
% espacio en coordenadas ECEF y en metros.

% A partir de los 100km de altura sobre la elipsoide devuelve siempre 1.

function [n] = refractiveIndex(x)

lla = ecef2llaGeod(x);
h = lla(3);

if(h < 0)
        %fprintf("La altura sobre la elipsoide es negativa. El indice de refraccion no puede ser calculado.\n");
        % n = 1;
         N = 77.604 * 1013.25/288.16;
         n = 1 + ( N / 1e6 );
        return
end

rho0 = 1013.25; % Presi�n atmosferica estandar a nivel del mar. [hPa]
g = 9.80665; % Aceleraci�n gravitacional en la superficie terrestre. [m/s^2]
T0 = 288.16; % Temperatura estandar a nivel del mar. [K]
M = 0.02896968; % Masa molar de aire seco. [kg/mol]
R0 = 8.314462618; % Constante de gas universal. [J/(mol�K)]

pressure = rho0 * exp (- (g * h * M)/(T0 * R0) );

    if (h < 11000) % Por debajo de los 11km de altura.
        temp = T0 - 6.5 * (h / 1000);
    elseif (h < 20000) % Por debajo de los 20km
        temp = 216.65; %[K]
    elseif (h < 32000)
        temp = 216.65 + 0.001 * (h - 20000);
    elseif (h < 48500)
        temp = 228.65 + 0.0025 * (h - 32000);
    elseif (h < 52500)
        temp = 270;
    elseif (h < 60500)
        temp = 270 - 0.00175 * (h - 52500);
    elseif (h < 80000)
        temp = 256 - (11/3250) * (h - 60500);
    elseif (h < 90500)
        temp = 190;
    elseif (h < 100000)
        temp = 190 + (3/2375) * (h - 90500);
    else
        N = 0; % A partir de esta altura, si bien la temperatura cambia la presi�n ya bajo tanto que se puede asumir N = 0
    end
        
    if (h < 100000)
        N = 77.604 * pressure/temp;
    end
    
    n = 1 + ( N / 1e6 );

    
    
end