%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% Script para correr los eventos desde el simulador %%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close 
% Se agrega a la carpeta el path donde estan los .csv del simulador
addpath("C:\Users\USUARIO\OneDrive - Alumnos Facultad de Ingeniería - UNLP\Tesina de Grado\GeneradorGPS\cmake-build-debug\")
addpath("C:\Users\USUARIO\OneDrive - Alumnos Facultad de Ingeniería - UNLP\\BECA Senyt\Repositorio Git Senyt\code-senyt\");
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trayectoria LEO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tierra_color();
POS_LEO= csvread('estimated_positions.csv');
POS_LEO =POS_LEO(:,2:4); 
hold on;    
for i = 10000:100:20000

    if i==10000
        plot3(POS_LEO(i,1),POS_LEO(i,2),POS_LEO(i,3), '+', 'MarkerSize', 10, 'color', 'magenta','LineWidth',2)
    end

    plot3(POS_LEO(i,1),POS_LEO(i,2),POS_LEO(i,3), '+', 'MarkerSize', 4, 'color', 'yellow','LineWidth',2)

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Trayectoria GPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
POS_GPS = csvread('propGPS.csv'); % PRN, TIEMPO, POS
POS_PRN = zeros(length(POS_GPS)/30,3);
PRN = 25; % SATÉLITES DISPONIBLES DEL 2 AL 32 SIN EL 28
cnt =1;
for i=1:length(POS_GPS)

    satelite = POS_GPS(i,1);

    if (satelite == PRN )

        POS_PRN(cnt,:) = POS_GPS(i,3:5);
        cnt = cnt + 1 ;
    end

end

for i = 1:length(POS_PRN)

    plot3(POS_PRN(i,1),POS_PRN(i,2),POS_PRN(i,3), '+', 'MarkerSize', 4, 'color', 'b','LineWidth',2)

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVENTOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grafica los eventos totales, no tiene que ser con el PRN que graficamos
% arriba.
% Para saber con que GPS se dio el evento podemos ver la tabla matriz
% eventos (matrizeventos).

hold on

P_GPS = csvread('POS_GPS_RO.csv');
P_LEO = csvread('POS_LEO_RO.csv');

for i = 1:length(P_GPS)

    plot3([P_LEO(i,1) P_GPS(i,1)],[P_LEO(i,2) P_GPS(i,2)],[P_LEO(i,3) P_GPS(i,3)], '-', 'MarkerSize', 8, 'MarkerFaceColor', 'r','LineWidth',2)

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CARACTERIZACIÓN DE EVENTOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
matriz_datos = csvread('duracion.csv');
vel_USAT = csvread('velocidad_LEO.csv');
close all
tierra_color
LEO_cont = matriz_datos(:,2:4);
GPS_cont = matriz_datos(:,5:7);
hold on

for i = 1:length(LEO_cont)/100

    %plot3([LEO_cont(i,1) GPS_cont(i,1)],[LEO_cont(i,2) GPS_cont(i,2)],[LEO_cont(i,3) GPS_cont(i,3)], '-', 'MarkerSize', 8, 'MarkerFaceColor', 'r','LineWidth',2)
    plot3(LEO_cont(i,1), LEO_cont(i,2), LEO_cont(i,3), '+', 'MarkerSize', 4, 'color', 'b','LineWidth',2)
end

% Reacomodamos por satélite
sort_events_LEO = zeros(size(LEO_cont));
sort_events_GPS = zeros(size(GPS_cont));
sort_vel_LEO = zeros(size(vel_USAT));
cnt = 1;
PRN_S = zeros(1,length(LEO_cont))'; 
for PRN = 1:33
    for i = 1:length(LEO_cont)
        if (matriz_datos(i,1) == PRN)

            sort_events_LEO(cnt,:) = matriz_datos(i,2:4);
            sort_events_GPS(cnt,:) = matriz_datos(i,5:7);
            sort_vel_LEO(cnt,:) = vel_USAT(i,:);
            PRN_S(cnt) = matriz_datos(i,1);
            cnt = cnt + 1;

        end
    end
end

% Calculo de duración de eventos, la matriz duracion_events tiene cuanto
% durael evento y a que satélite le corresponde  y en que satélite 
% El archivo p interes te dice cuantos evetnos tiene que haber 

step = norm(sort_events_LEO(1,:)-sort_events_LEO(2,:))*1.2;
cnt = 1;
duracion = 1;
for i = 1:length(LEO_cont)-1



    if norm(sort_events_LEO(i,:)-sort_events_LEO(i+1,:)) < step
    
        duracion = duracion + 1;
    
    else 

    duracion_events(cnt,1) = duracion;
    duracion_events(cnt,2) = PRN_S(i);
    cnt = cnt + 1;
    duracion = 1;

    end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANTES DE PROPAGAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Necesitamos tener una representación de cual es el evento de RO para cada
% par, asi tenemos una posición.

% Tomamos como criterio elegir el que tenga una altitud menor antes de
% cambiar de signo. 
step = norm(sort_events_LEO(1,:)-sort_events_LEO(2,:))*1.2;

h_comp = 10e10; % Altura de comparación
colum = 1;
cnt = 1;
clear altitudes
for i = 1:length(LEO_cont)-1

    if (norm(sort_events_LEO(i,:)-sort_events_LEO(i+1,:))<step)
        % Cálculo de la altura minima entre el rayo GPS-USAT
        l = sort_events_GPS(i,:) - sort_events_LEO(i,:);
        p = -dot(l,sort_events_LEO(i,:))/norm(l)^2;
        minH = p*l+sort_events_LEO(i,:);
        minHlla = ecef2llaGeod(minH);
        hminLOS = minHlla(3);

        % Listado de altitudes minimas por satélite y la posicion del
        % satélite en forma de índice

        altitudes(cnt,[colum, colum+1]) = [hminLOS i];

        cnt = cnt + 1;

    else
        colum = colum + 2;
        cnt = 1;

    end



end

% Altura contiene una fila de altitudes y una fila del indice que tengo que
% ir a buscar a sort_event para tener la posicion del LEO y el GPS que
% tienen la menor altitud. 

% Encuentro la altitud más chica en cada evento 
tierra_color
cnt = 1;
[filas,columnas] = size(altitudes);
for i = 1:columnas/2

    impares = 2*i - 1;

    % Saco de la matriz los 0 y pongo todo en valor absoluto
    filteredData = abs(altitudes(abs(altitudes(:,impares)) ~= 0, impares));
    % Obtengo la altitud mínima y el índice que tengo que ir a buscar
    % mi posición de USAT y GPS.
    [alt, minIndex] = min(filteredData);
    
    minIndex = altitudes(minIndex,impares + 1); % Me llevo el indice correspondiente

    %
    RO_events_USAT(:,cnt) = sort_events_LEO(minIndex,:)';
    RO_events_GPS(:,cnt)  = sort_events_GPS(minIndex,:)';
    RO_events_USAT_vel(:,cnt) = sort_vel_LEO(minIndex,:)';
    hold on
    plot3([RO_events_USAT(1,cnt) RO_events_GPS(1,cnt)],[ RO_events_USAT(2,cnt)  RO_events_GPS(2,cnt)],[ RO_events_USAT(3,cnt) RO_events_GPS(3,cnt)], '-', 'MarkerSize', 8, 'MarkerFaceColor', 'r','LineWidth',2)
    cnt = cnt + 1;

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ÁNGULOS Y EJES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Si corremos solo esto caargar posiciones de satélites
% EJES DEFINIDOS SOBRE EL SATÉLITE
hold on

% close all;

tierra_color()

for i=10:10

    % Con Event determinamos cual de los eventos de los 403 queremo ver
    Event = i;
    plot3([RO_events_USAT(1,Event) RO_events_GPS(1,Event)],[ RO_events_USAT(2,Event)  RO_events_GPS(2,Event)],[ RO_events_USAT(3,Event) RO_events_GPS(3,Event)], 'o', 'MarkerSize', 18, 'MarkerFaceColor', 'r','LineWidth',2)

    % Eje desde el centro de la tierra al LEO

    a1 = -1*RO_events_USAT(:,Event);

    a1 = a1/norm(a1); % Lo hacemos unitario

    % Eje perpendicular determinado por la velocidad

    a2 = RO_events_USAT_vel(:,Event);
    a2 = a2/norm(a2);

    % Eje perpendicular a la velocidad y la posición

    a3 = cross(a1,a2);
    a3 = a3/norm(a3);

    % Dibujamos los ejes en el origen
    % quiver3(0,0,0,a1(1)*10e6,a1(2)*10e6,a1(3)*10e6, 'LineWidth', 2,'Color','blue') % eje x
    % quiver3(0,0,0,a2(1)*10e6,a2(2)*10e6,a2(3)*10e6, 'LineWidth', 2,'Color','red') % eje y
    % quiver3(0,0,0,a3(1)*10e6,a3(2)*10e6,a3(3)*10e6, 'LineWidth', 2,'Color','black') % eje z

    % Dibujamos los ejes en el USAT
    quiver3(RO_events_USAT(1,Event),RO_events_USAT(2,Event),RO_events_USAT(3,Event),a1(1)*10e6,a1(2)*10e6,a1(3)*10e6, 'LineWidth', 2,'Color','blue') % eje x
    quiver3(RO_events_USAT(1,Event),RO_events_USAT(2,Event),RO_events_USAT(3,Event),a2(1)*10e6,a2(2)*10e6,a2(3)*10e6, 'LineWidth', 2,'Color','red') % eje y
    quiver3(RO_events_USAT(1,Event),RO_events_USAT(2,Event),RO_events_USAT(3,Event),a3(1)*10e6,a3(2)*10e6,a3(3)*10e6, 'LineWidth', 2,'Color','black') % eje z

    % Obtenemos la dirección final del disparo al resolver la ecuación
    % diferencial que modela la señal
    % quiver3(0,0,0,y(1,end)*10e6,y(2,end)*10e6,y(3,end)*10e6, 'LineWidth', 2,'Color','magenta') % eje z



    % Resuelvo la señal de RO
    [x,y,y0] = propagador_de_senial_RO(i);
    % graficar_x(x');
    % Cálculo de ángulos de elevación y azimuth
    quiver3(RO_events_USAT(1,Event),RO_events_USAT(2,Event),RO_events_USAT(3,Event),-1*y(1,end)*10e6,-1*y(2,end)*10e6,-1*y(3,end)*10e6, 'LineWidth', 2,'Color','magenta') % eje z

    % Ángulo de inclinación en primer eje (TIENE QUE ESTAR ENTRE +-55)
    [theta,flag] = angulos_para_antenna(a2*10e6,-1*y(:,end)*10e6,a2*10e6,a1*10e6);
    if flag == 1
        angulos(i,1) = -theta; % puede ser que dé negativo...
    else
        angulos(i,1) = theta;
    end
    
    % Ángulo de inclinación en segundo eje (TIENE QUE ESTAR ENTRE MAS MENOS 20)
    [theta,flag] = angulos_para_antenna(a2*10e6,-1*y(:,end)*10e6,a2*10e6,a3*10e6);
    if flag == 1
        angulos(i,2) = -theta;
    else
        angulos(i,2) = theta;
    end


end
%% Posición óptima de antena

% grid on
% hold on;plot(angulos(:,1),angulos(:,2),'*')
% xlabel('Elv','Interpreter','latex');
% ylabel('Az','Interpreter','latex');

lim_inf_inc = 0;
lim_sup_inc = 50;

lim_inf_az = 0;
lim_sup_az = 20;

contador = 1;

fila = 1;

while(lim_sup_inc<90 || lim_sup_az< 90)

lim_inf_inc = lim_inf_inc +1;
lim_sup_inc = lim_sup_inc +1;

lim_inf_az = lim_inf_az +1;
lim_sup_az = lim_sup_az +1;
   
for i=1:length(angulos)
    
    if (angulos(i,1)< lim_sup_inc && angulos(i,1)> lim_inf_inc && angulos(i,2)< lim_sup_az && angulos(i,2)> lim_inf_az)

    contador = contador +1;
    end
end

matrix(fila,[ 1 2 3 4 5]) = [lim_inf_inc lim_sup_inc lim_inf_az lim_sup_az contador ];

fila = fila +1;
contador =1;

end

%% Simulación / Demo

% Correr el primer bloque y cargar las pos y vel


% La idea es hacer la orbita completa ver el evento, propagar, sacar la
% dirección plotear los ejes. 

% Tomamos la órbita del USAT-1       

tierra_color();
POS_LEO= csvread('estimated_positions.csv');
POS_LEO =POS_LEO(:,2:4); 
hold on;    

for i = 1:50:2600

    if i==10000
        plot3(POS_LEO(i,1),POS_LEO(i,2),POS_LEO(i,3), '-o', 'MarkerSize', 10, 'color', 'magenta','LineWidth',2)
    end

    plot3(POS_LEO(i,1),POS_LEO(i,2),POS_LEO(i,3), '-o', 'MarkerSize', 4, 'color', 'yellow','LineWidth',2)
    
    pause(0.05)

end

% Detectamos el evento (GRAFICAMOS AMBAS POSICIONES USAT-GPS)

plot3(RO_events_USAT(1,1),RO_events_USAT(2,1),RO_events_USAT(3,1), '-o', 'MarkerSize', 7, 'color', 'RED','LineWidth',2)
plot3(RO_events_GPS(1,1),RO_events_GPS(2,1),RO_events_GPS(3,1), '-o', 'MarkerSize', 7, 'color', 'RED','LineWidth',2)

pause(2)
plot3([RO_events_USAT(1,1) RO_events_GPS(1,1)],[RO_events_USAT(2,1) RO_events_GPS(2,1)],[RO_events_USAT(3,1) RO_events_GPS(3,1)], '-', 'MarkerSize', 4, 'MarkerFaceColor', 'r','LineWidth',2)
pause(2)
% Resolvemos la trayectoria 
[x,y,y0] = propagador_de_senial_RO(1);

graficar_x(x')
pause(2)
% Se ponen los ejes y la dirección final
Event=1;
a1 = -1*RO_events_USAT(:,Event);

a1 = a1/norm(a1); % Lo hacemos unitario 

% Eje perpendicular determinado por la velocidad 

a2 = RO_events_USAT_vel(:,Event);
a2 = a2/norm(a2);

% Eje perpendicular a la velocidad y la posición 

a3 = cross(a1,a2);
a3 = a3/norm(a3);
quiver3(RO_events_USAT(1,Event),RO_events_USAT(2,Event),RO_events_USAT(3,Event),a1(1)*10e6,a1(2)*10e6,a1(3)*10e6, 'LineWidth', 2,'Color','blue') % eje x
quiver3(RO_events_USAT(1,Event),RO_events_USAT(2,Event),RO_events_USAT(3,Event),a2(1)*10e6,a2(2)*10e6,a2(3)*10e6, 'LineWidth', 2,'Color','red') % eje y 
quiver3(RO_events_USAT(1,Event),RO_events_USAT(2,Event),RO_events_USAT(3,Event),a3(1)*10e6,a3(2)*10e6,a3(3)*10e6, 'LineWidth', 2,'Color','black') % eje z

pause(2)
quiver3(RO_events_USAT(1,Event),RO_events_USAT(2,Event),RO_events_USAT(3,Event),-1*y(1,end)*10e6,-1*y(2,end)*10e6,-1*y(3,end)*10e6, 'LineWidth', 2,'Color','magenta','MaxHeadSize',3) % eje z