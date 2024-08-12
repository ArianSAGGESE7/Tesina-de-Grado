%% Propagación de orbita USAT

% Un eje quiero que sea apuntando hacia la tierra
% El otro para un recorrido con exentricidad nula es la dirección
% tangencial
% Con esos dos me queda definida la terna


% Generamos algunas orbitas de usat 
clc;clear all;close all


G = 6.67430e-11; % Constante de gravitación universal (m^3/kg/s^2)
M = 5.972e24;    % Masa de la Tierra (kg)
R = 6366000;     % Radio de la Tierra (m)
 
dataLEO = yumaread("usat.yuma.txt");
tLEO = dataLEO.Time(1);

Periodo=2*pi*sqrt((dataLEO.SQRTA^2)^3/(M*G));
resolucion = 10 ; % [segundos]
vueltas = 1;


for i=1:(Periodo/resolucion)*vueltas+resolucion

    [satPos,satVel] = gnssconstellation(tLEO,dataLEO,GNSSFileType="YUMA");

    reg_tLEO(i) = tLEO;

    tLEO.Second=tLEO.Second+resolucion;

    Posicion_LEO(:,i) = satPos;
    Velocidad_LEO(:,i) = satVel;

end

% tierra_color()
% comet3(Posicion_LEO(1, :), Posicion_LEO(2, :), Posicion_LEO(3, :))

 %% Ejes de USAT con graficos  (hecho con la última posición USAT) (saque los ejes con la velcidad {tangencial})

close all
tierra_color()


for p=1:10:length(Posicion_LEO)
pos =p;
% Vector desde el LEO hasta el centro de la tierra
a1 = Posicion_LEO(:,pos);
b1 = zeros(1,3)';

% Calcula el vector dirección del vector que va de a a b
dir_vector1 = (b1 - a1)/2;

% Dibuja el vector quiver3
hold on;
quiver3(a1(1), a1(2), a1(3), dir_vector1(1), dir_vector1(2), dir_vector1(3), 'LineWidth', 2,'Color','blue');
plot3(a1(1), a1(2), a1(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Punto a
plot3(b1(1), b1(2), b1(3), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Punto b

% Vector tangente a la trayectoria que me termina de definir la terna

a2 = Posicion_LEO(:,pos);

b2 = Velocidad_LEO(:,pos)+Posicion_LEO(:,pos);

% Calcula el vector dirección del vector que va de a a b
dir_vector2 = (b2 - a2)*500;

quiver3(a2(1), a2(2), a2(3), dir_vector2(1), dir_vector2(2), dir_vector2(3), 'LineWidth', 2,'color','red');

% Vector perpendicular a Dirvector y Dirvector1

c = cross(dir_vector1,dir_vector2);
c=c/norm(c);
a3 = Posicion_LEO(:,pos);
b3 = c + Posicion_LEO(:,pos);
dir_vector3 = (b3 - a3)*3e6;


quiver3(a3(1), a3(2), a3(3), dir_vector3(1), dir_vector3(2), dir_vector3(3), 'LineWidth', 2,'color','black');

% Patrón simbólico de radiación de la antena

P = -1*dir_vector1;  % Center point of the circle
V = dir_vector1/10;  % Vector defining the normal plane
R = tan(30*pi/180)*norm(-1*dir_vector1-a1);          % Ingresamos el ángulo
hold on;
[x_circle, y_circle, z_circle] = plot_cone(P, V, R);
quiver3(P(1), P(2), P(3), V(1), V(2), V(3), 'g');  % Mostrar el vector V
plot3(x_circle, y_circle, z_circle, 'b-');
scatter3(-1*dir_vector1(1), -1*dir_vector1(2), -1*dir_vector1(3), 'magenta', 'filled');
scatter3(x_circle(1), y_circle(1),z_circle(1), 'magenta', 'filled');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Circle in the Normal Plane Defined by Vector V');
grid on;

P = a1;

for i = 1:length(x_circle)
    plot3([P(1), x_circle(i)], [P(2), y_circle(i)], [P(3), z_circle(i)], 'b-');
end

pause(1)
end