% -------------------------------------------------------
% Modelado de la señal de RO mediante pantallas de fase -
% -------------------------------------------------------

% Idea: modelar la atmósfera como pantallas de fase que cada una modifica el número de onda
% Para ello hay que definir algunos parámetros de interés

% Parámetros
close all; clc;clear all
RE = 6371e3;
lambda = 0.19029; % Longitud de onda de la señal L1 GPS 
k = 2 * pi / lambda; % Número de onda
num_screens = 2000; % Número de pantallas de fase
delta_l = 1e3; % Separación entre pantallas 
H = 524.288e3; % Extensión vertical de la pantalla
M = 2^10;
z = linspace(H/2, -H/2,M); % Puntos de muestreo vertical
N = @(h) 400 * exp(-h/8e3); % Refractividad modelo A (exponencial)

% Para cada una de las pantallas calculamos el exceso de fase

for i=1:num_screens
    
    l = num_screens/2*delta_l - delta_l*(i-1); % Acá ponemos la distancia en "y" que nos movemos
    cos_theta = (1-l^2/RE^2)^(1/2);
    z_sh = RE*(cos_theta - 1); % Altura
    h = (z-z_sh)/cos_theta; % Altura en la que vamos a valuar 
    s(:,i) = 1e-6*N(h)*delta_l; % Exceso de fase por cada delta_l en todo z 

end

% Ahora que teneos el exceso de fase para cada "l" en todo los "z" podemos
% ir a la formula cerrada de la solución de la ecuación y calcular el
% aporte de todos los rayos por medio de la FFT

% Parámetros del modelo de atenuación
alpha_0 = 10e-4; % Coeficiente de atenuación inicial 
H_0 = 8e3; % Escala de altura 

% Modelo de atenuación atmosférica
Windows_atenuacion = exp(-alpha_0 * exp(-z/ H_0)); % hay que calibrarla para que quede acorde al paper

% Inicializa la señal de entrada (asumiendo una señal plana de entrada)
U_in = ones(size(z));
%  a través de las pantallas
y_s = -num_screens/2*delta_l;
for i = 1:num_screens

    % Calcula la amplitud compleja de salida en la pantalla actual (un instante luego)

    U_out = U_in .* exp(1i * k * s(:, i)).* Windows_atenuacion';
    
    % Fourier (acá hago la suma tal cual el paper)
    U_out_fft = fft(U_out);
    
    % Calcula la señal en el plano de observación (propagación libre) (tengo que ponerle la parte de
    % corrimiento en los planos en y)
    ky = sqrt(k^2 - (2 * pi *  (0:length(z) - 1)  / H).^2);
    u_observed(:,i) = ifft(U_out_fft(:,1) .* exp(1i* ky' *delta_l));
    % Actualizao para la sig pantalla
    U_in = u_observed;
    
end



% Actualiza la señal de entrada para la siguiente pantalla
U_in = u_observed;

% Amplitud y fase de la señal observada en el LEO
amplitud_LEO = abs(u_observed);
phase_LEO = angle(u_observed);

close all
% Gráficas de la señal
figure;
%subplot(2, 1, 1);
plot(z/1000, amplitud_LEO(:,end)');
xlabel('Altura (km)');
ylabel('Amplitud');
title('Amplitud de la señal en el LEO');
%%
subplot(2, 1, 2);
plot(z/1000, phase_LEO);
xlabel('Altura (km)');
ylabel('Fase (radianes)');
title('Fase de la señal RO en LEO');
