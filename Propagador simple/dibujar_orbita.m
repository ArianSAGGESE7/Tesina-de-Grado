% Este codigo es para entender la mecanica orbital bajo un modelo simple de
% dos cuerpos 
% Considerando solo el efecto de la tierra sobre el satélite

% Definir el radio de la Tierra según el modelo WGS84

radio_tierra = 6371e3; % en kilómetros
escala_plano = 2;


 % [x, y, z] = sphere(100); 
 % earth_surface = surf(x * radio_tierra, y * radio_tierra, z * radio_tierra);
tierra_color();
hold on;

theta = linspace(0, 2 * pi, 100);
x_ecuador = radio_tierra * cos(theta);
y_ecuador = radio_tierra * sin(theta);
z_ecuador = zeros(size(theta)); % el ecuador está en el plano z=0

ecuador_coords = [x_ecuador' * escala_plano, y_ecuador' * escala_plano, z_ecuador'];

% Trazar el plano del ecuador extendido
fill3(ecuador_coords(:,1), ecuador_coords(:,2), ecuador_coords(:,3), 'r');

% Configuración del aspecto de la figura
axis equal; % Para que los ejes tengan la misma escala
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
quiver3(0, 0, 0, 2*radio_tierra, 0, 0, 'k', 'LineWidth', 2); % Eje x
quiver3(0, 0, 0, 0, 2*radio_tierra, 0, 'g', 'LineWidth', 2); % Eje y
quiver3(0, 0, 0, 0, 0, 2*radio_tierra, 'b', 'LineWidth', 2); % Eje z


% Definir parámetros orbitales
a = 6371e3+600e3; % Semieje mayor en kilómetros
e = 0; % Excentricidad
i = deg2rad(97.77); % Inclinación en radianes
Omega = deg2rad(0); % Longitud del nodo ascendente en radianes
omega = deg2rad(0); % Argumento del perigeo

% Crear un vector de anom alías verdaderas
nu = linspace(0, 2*pi, 1000);

% Calcular las coordenadas orbitales en el sistema orbital
r_orbital = a * (1 - e^2) ./ (1 + e * cos(nu));

x_p = r_orbital.*cos(nu + omega);
y_p = r_orbital.*sin(nu + omega);
x_s = x_p.*cos(Omega)-y_p.*cos(i).*sin(Omega);
y_s = x_p.*sin(Omega)+y_p.*cos(i).*cos(Omega);
z_s = y_p.*sin(i);

% Plotear la órbita
hold on;
plot3(x_s, y_s, z_s, 'b', 'LineWidth', 2);


axis equal; % Para que los ejes tengan la misma escala
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Órbita alrededor de la Tierra');

% Ajustar la vista para que se muestre correctamente
grid on;

[x,v] = kepler2eci(a, e, i, Omega, omega, pi/3);
plot3(x(1),x(2),x(3),'MarkerEdgeColor',"r",'LineWidth',10,'MarkerSize',6,'Marker','*')


legend('Orbitad del Usat-1','interpreter')
