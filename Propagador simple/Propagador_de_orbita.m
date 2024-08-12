function [x,v]=Propagador_de_orbita(a, e, i, Omega, omega, nu,Periodo,resolucion,vueltas)

% FUNCION PARA ESTIMAR LA ORBITA DADOS LOS PARÁMETROS KEPLERIANOS 

%ENTRADAS

% a semieje mayor
% e excentricidad [°]
% Omega Longitud del nodo ascendente [°]
% omega Argumento del perigeo [°]
% i inclinación [°]
% nu Anomalia verdadera [°]
% Periodo = 2*pi*sqrt(a^3/(mu))  % MU =3.985891960000000e+14
% Resolución =  paso de integración gralmente 1s
% vueltas = cantidad de vueltas gralmente 1



G = 6.67430e-11; % Constante de gravitación universal (m^3/kg/s^2)
M = 5.972e24;    % Masa de la Tierra (kg)
R = 6371000;     % Radio de la Tierra (m)
dt = resolucion;          % Paso de tiempo (s)
num_steps = floor(Periodo/dt)*vueltas+dt; % Número de pasos de tiempo

% Definimos posición y velocidad

% [x0] = kepler2ecef(R+1500e3, 0, 45, 30, 30, 20); % Parámetros de órbita

[x0,v0] = kepler2eci(a, e, i, Omega, omega, nu); % Parámetros de órbita
% [x1] = kepler2ecef(a, e, i, Omega, omega, nu+0.00001);% Utilizamos un punto muy cercano para calcular la direccion de la velocidad

% dir_tangente = x0-x1;
% dir_tangente=-1*dir_tangente/(norm(dir_tangente)).*sqrt(G*M/(norm(x0)));
% 
% v0 = dir_tangente;


aceleracion = @(x, y, z) -G * M / ((x^2 + y^2 + z^2)^(3/2));

grad = @(t, state) [state(4); state(5); state(6); ...
    aceleracion(state(1), state(2), state(3)) * state(1); ...
    aceleracion(state(1), state(2), state(3)) * state(2); ...
    aceleracion(state(1), state(2), state(3)) * state(3)];

time = zeros(num_steps, 1);
states = zeros(num_steps, 6);

state = [x0; v0];

% Runge-kutta de orden 4 para la resolución de la ecuación diferencial
for i = 1:num_steps

    time(i) = (i - 1) * dt;
    states(i, :) = state';

    k1 = grad(time(i), state);
    k2 = grad(time(i) + dt/2, state + (dt/2) * k1);
    k3 = grad(time(i) + dt/2, state + (dt/2) * k2);
    k4 = grad(time(i) + dt, state + dt * k3);

    % Actualizar el estado utilizando RK4
    state = state + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% Graficar la trayectoria tridimensional del satélite
 
tierra_color();
hold on
comet3(states(:, 1), states(:, 2), states(:, 3))
x = states(:,1:3);
v = states(:,4:6);
end

