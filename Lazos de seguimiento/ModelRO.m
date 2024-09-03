% Modelado de la señal de RO mediante pantallas de fase

% Idea modelar la atmósfera como pantallas de fase que cada una modifica el número de onda
% Para ello hay que definir algunos parámetros de interés

% División entre frentes (según el paper se realizan 200)

% Parámetros
close all; clc
R = 6371e3;
lambda = 0.19029; % Longitud de onda de la señal L1 GPS (~19 cm)
k = 2 * pi / lambda; % Número de onda
num_screens = 2000; % Número de pantallas de fase
delta_l = 1e3; % Separación entre pantallas (1 km)
H = 524.288e3; % Extensión vertical de la pantalla (524.288 km)
M = 2^10;
z_values = linspace(H/2, -H/2,M); % Puntos de muestreo vertical

l_0 = 750e3+R; % Separación desde el centro de la tierra hasta el satélite LEO 
% Perfil de refractividad (ejemplo simplificado)
N = @(h) 400 * exp(-h/8e3); % Refractividad modelo A (exponencial)

% Calcula la fase excedente para cada pantalla

% Tener en cuenta que el problema se plantea de manera bi-dimensional
% Por lo que "y" es el avance de la onda y "z" representa el movimiento en vertical

excess_phase = zeros(length(z_values), num_screens);
for i = 1:num_screens
    h = z_values + (i - 1) * delta_l; % Altura correspondiente a la pantalla
    refractivity = N(h); % Refractividad en esa altura
    excess_phase(:, i) = 1e-6 * refractivity * delta_l; % Fase excedente
end

% Parámetros del modelo de atenuación
alpha_0 = 1e-4; % Coeficiente de atenuación inicial (valor de ejemplo)
H_0 = 4e3; % Escala de altura (ejemplo para la troposfera húmeda)

% Modelo de atenuación atmosférica
attenuation_factor = exp(-alpha_0 * exp(-z_values / H_0));


% Inicializa la señal de entrada (asumiendo una señal plana de entrada)
U_in = ones(size(z_values));

% Propaga la señal a través de las pantallas
for i = 1:num_screens
    % Calcula la amplitud compleja de salida en la pantalla actual (un instante luego)
    U_out = U_in .* exp(1i * k * excess_phase(:, i)).* attenuation_factor';
    
    % Realiza la expansión en Fourier (acá hago la suma tal cual el paper)
    U_out_fft = fft(U_out);
    
    % Calcula la señal en el plano de observación (propagación libre) (tengo que ponerle la parte de
    % corrimiento en los planos en y)
    ky =sqrt(k^2 - (2 * pi * (0:length(z_values) - 1) / H).^2);
    u_observed = ifft(U_out_fft .* exp(1i * ky' * delta_l));
    
    % Actualiza la señal de entrada para la siguiente pantalla
    U_in = u_observed;
end


% Actualiza la señal de entrada para la siguiente pantalla
U_in = u_observed;


% Amplitud y fase de la señal observada en el LEO
amplitude_LEO = abs(u_observed);
phase_LEO = angle(u_observed);


close all
% Gráficas de la señal
figure;
subplot(2, 1, 1);
plot(z_values/1000, amplitude_LEO/max(amplitude_LEO));
xlabel('Altura (km)');
ylabel('Amplitud');
title('Amplitud de la señal en el LEO');

subplot(2, 1, 2);
plot(z_values/1000, phase_LEO);
xlabel('Altura (km)');
ylabel('Fase (radianes)');
title('Fase de la señal RO en LEO');

% Calcula la transformada de Fourier de la señal final observada en el LEO
U_fft = fftshift(fft(u_observed));

fs = 1 / delta_l; % Frecuencia de muestreo en el espacio
f = linspace(-fs/2, fs/2, length(U_fft));

U_ftt_abs = abs(U_fft);

% Gráfico del espectro de la señal en el LEO
figure;
plot(f, U_ftt_abs);
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
title('Espectro de la señal RO en LEO');
grid on;


% Cálculo de aceleración de fase
delta_t = 1/50; % Lo tomamos del paper
for i = 2:length(phase_LEO)-1
    diff_phase(i) = (phase_LEO(i+1)-2*phase_LEO(i)+phase_LEO(i-1))/(2*pi*delta_t^2);
end
figure;
plot(z_values(2:end),diff_phase);
title("Aceleración de fase")
