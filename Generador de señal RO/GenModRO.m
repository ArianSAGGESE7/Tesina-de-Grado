% ================================================
%  Perfiles de fase y atenuación para señal de RO
% ================================================

clear all;clc;close all
% Parámetros
Re = 6371e3; % Radio terrestre 
lambda = 0.19; % Longitud de onda de la señal L1 GPS (~19 cm)
k = 2 * pi / lambda; % Número de onda
num_screens = 2000; % Número de pantallas de fase
delta_l = 1e3; % Separación entre pantallas (1 km)
H = 524.288e3; % Extensión vertical de la pantalla (524.288 km)
M = 2^10;
z_values = linspace(-H/2, H/2, M); % Puntos de muestreo vertical
cosTheta = @(l) (1-l^2/Re^2)^(0.5); % Angulo de proyección en el mapa 

% Perfil de refractividad (ejemplo simplificado)
N = @(h) 400 * exp(-h/8e3); % Refractividad modelo A (exponencial)

% Parámetros del modelo de atenuación
alpha_0 = 3e-4; % Coeficiente de atenuación inicial 
H_0 = 60e3; % Escala de altura      

% Modelo de atenuación atmosférica
a_w= exp(-alpha_0 * exp(-z_values/ H_0))';

% Calcula la fase excedente para cada pantalla
excess_phase = zeros(length(z_values), num_screens);
for i = 1:num_screens

    l =  -num_screens/2*delta_l + delta_l*i; % me situo en el y_S
    z_sh = Re*(cosTheta(l)-1); % calculo el z_sh
    h = (z_values-z_sh)./cosTheta(l); % Altura correspondiente a la pantalla
    refractivity = N(h); % Refractividad en esa altura
    excess_phase(:, i) = 1e-6 * refractivity * delta_l; % Fase excedente

    % Parámetros del modelo de atenuación
    alpha_0 = 3e-4; % Coeficiente de atenuación inicial 
    H_0 = 40e3; % Escala de altura 

    % Modelo de atenuación atmosférica
    a_w(:,i)= exp(-alpha_0 * exp(-h/ H_0))';

end


% Inicializa la señal de entrada (asumiendo una señal plana de entrada)
U_in = ones(length(z_values),1);

% Propaga la señal a través de las pantallas
for i = 1:num_screens
    % Calcula la amplitud compleja de salida en la pantalla actual
    U_out = a_w(:,i).*U_in .* exp(1i * k * excess_phase(:, i));
    
    % Realiza la expansión en Fourier
    U_out_fft = fft(U_out);
    % Calcula la señal en el plano de observación (propagación libre)
    ky = sqrt(k^2 - (2 * pi * (-M/2:M/2-1) / H).^2);
    u_observed = ifft(U_out_fft .* exp(1i * ky' * (3000e3-delta_l*i)));
    
    % Actualiza la señal de entrada para la siguiente pantalla
    U_in = u_observed;
end



% Amplitud y fase de la señal observada en el LEO
amplitude_LEO = abs(u_observed);
phase_LEO = angle(u_observed);

% Gráficas de la señal
figure;
subplot(2, 1, 1);
plot(z_values, amplitude_LEO/max(amplitude_LEO));
xlabel('Altura (m)');
ylabel('Amplitud');
title('Amplitud de la señal RO en LEO');

subplot(2, 1, 2);
plot(z_values, phase_LEO);
xlabel('Altura (m)');
ylabel('Fase (radianes)');
title('Fase de la señal RO en LEO');


% Cálculo de aceleración de fase
delta_t = 1/50; % Lo tomamos del paper
for i = 2:length(phase_LEO)-1
    diff_phase(i) = (phase_LEO(i+1)-2*phase_LEO(i)+phase_LEO(i-1))/(2*pi*delta_t^2);
end
figure;
plot(z_values(2:end),diff_phase);
title("Aceleración de fase")

% Doppler 

doppler = diff(phase_LEO);
figure;
plot(z_values(2:end),diff_phase);
title("Doppler")