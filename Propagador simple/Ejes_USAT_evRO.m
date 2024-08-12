%% Script de ejes con disparos propios de RO
% La idea de este script es cálcular el ángulo del punto final de la
% resolución de cada uno de los eventos de RO que simulamos, con respecto a
% una terna espacial centrada en el USAT que apunta a la tierra y con el
% patrón de ganancia de antena como eje x+. Con esto en mente...

% Cargamos las posiciones de USAT Y GPS para poder hacer en cada una de ellas una simulación
% gráfica de los ejes y el patrón de radiación que si bien no es un cono se
% asemeja a ello (respetamos el ángulo de azimuth)

 % close all;clear all;clc
 tierra_color()
% Carga de variables de otras simulaciones

rGPS = load('POSICIONES_PRN.mat');
rLEO = load('POSICIONES_USAT.mat');
velLEO = load('velocidad_usat.mat');
Posicion_final_ray = load('ultima_posicion.mat');
Exitosos = load('exitosos.mat');


rGPS=rGPS.EVENTO_PRN_fino;
Posicion_LEO = rLEO.EVENTO_USAT_fino;
Velocidad_LEO = velLEO.EVENTO_PRN_vel;
Posicion_final_ray = Posicion_final_ray.ultima_posicion;
Exitosos=Exitosos.exitosos; % contiene cuales de las posiciones de USAT son las que se pudieron resolver con los scripts de RO

counting = 1;
for p=1
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
    scatter3(rGPS(1,p), rGPS(2,p), rGPS(3,p),'magenta', 'filled','LineWidth',4)
    
    % -------------------------------------------------
    % Patrón simbólico de radiación de la antena
    % -------------------------------------------------
    % P = -1*dir_vector2;  % Center point of the circle
    % V = dir_vector2/10;  % Vector defining the normal plane
    % R = tan(60*pi/180)*norm(-1*dir_vector2-a1);          % Ingresamos el ángulo
    % hold on;
    % [x_circle, y_circle, z_circle] = plot_cone(P, V, R);
    % quiver3(P(1), P(2), P(3), V(1), V(2), V(3), 'g');  % Mostrar el vector V (vector normal)
    % plot3(x_circle, y_circle, z_circle, 'b-');

    % scatter3(-1*dir_vector1(1), -1*dir_vector1(2), -1*dir_vector1(3),
    % 'magenta', 'filled'); %punto sobre el eje de maxima ganancia
    % scatter3(x_circle(1), y_circle(1),z_circle(1), 'magenta', 'filled');
    % %punto en la esfera
    % axis equal;
    % grid on;
    % 
    % P = a1;
    % scatter3(Posicion_final_ray(1,counting), Posicion_final_ray(2,counting),Posicion_final_ray(3,counting), 'magenta', 'filled'); % grafico el último punto de RO
    % % incremento el contador de la variable de posición final del rayo
    % counting = counting + 1;
    % for i = 1:length(x_circle)
    %     plot3([P(1), x_circle(i)], [P(2), y_circle(i)], [P(3), z_circle(i)], 'b-');
    % end

    % pause(3) % Cambiar la pausa entre evento y evento
end