% Código de simulación para la lectura de semanas de GPS simulando una
% semana completa.
clc;clear all; close all
%%                                1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulación de constelación GPS + LEO orbita USAT1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Algunas variables de uso frecuente

G = 6.67430e-11; % Constante de gravitación universal (m^3/kg/s^2)
M = 5.972e24;    % Masa de la Tierra (kg)
R = 6366000;     % Radio de la Tierra (m)

% Leemos los datos del archivo Yb UMA
data = yumaread("almanac.yuma.week0100.061440.txt");
% data = yumaread("usat.yuma.txt");
% Resolvemos la orbita para un tiempo de simulación, una resolución
% temporal y una cantidad de vueltas enteras, es decir períodos

% gnssconstellation es un afuncion de matlab que interpreta los parametros
% Yuma y devuelve una aproximación de la posición y velocidad del satélite
% en coordenadas ECEF (Posición y velocidad)

t = data.Time(1);
% [satPos,satVel,satID] = gnssconstellation(t,data,GNSSFileType="YUMA");

% Una vez que tenemos la posición de la simulación de un dia en una
% constelación GPS simulamos para distintos tiempos de manera tal de
% recorrer la orbita.

% La resolución temporal es importante ya que determina el costo
% computacional


% Hay que determinar cuanto tiempo se pretende simular para completar una
% órbita.

% Generalmente una órbita GPS tiene un periodo de 11.58 hs, por lo que si
% le damos 2 vueltas tendremos el registro de 1 día, en posición,
% comenzando en la fecha que contiene la variable "data.time"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametro de satélite a simular (luego deberian ser todos)
SVn=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probamos para un satélite en particular PRN 1

datos_PRN1 = data(SVn,:);

% Posición inicial

[satPos,satVel1,satID] = gnssconstellation(t,datos_PRN1,GNSSFileType="YUMA");

Periodo=2*pi*sqrt((datos_PRN1.SQRTA^2)^3/(M*G));
resolucion = 10 ; % [segundos]
vueltas = 2;
% % Tengo que variar dentro de la estructura de tiempo los segundos
% % si la resolución temporal la tengo en segundos sino deberia cambiar
% 
% t.Second=t.Second+resolucion;
% 
% for i=1:floor(Periodo/resolucion*vueltas)+resolucion
% 
%     [satPos] = gnssconstellation(t,datos_PRN1,GNSSFileType="YUMA");
% 
%     t.Second=t.Second+resolucion;
% 
%     reg_tGPS(i) = t;
% 
%     Posicion_PRN1(:,i) = satPos';
% end


%
tierra_color()
% comet3(Posicion_PRN1(1, :), Posicion_PRN1(2, :), Posicion_PRN1(3, :))
%
% Ground_track(Posicion_PRN1)% chequeo track earth
%%
% Leemos el archivo LEO
dataLEO = yumaread("usat.yuma.txt");
tLEO = dataLEO.Time(1);


for i=1:(Periodo/resolucion)*vueltas+resolucion

    [satPos,satVel] = gnssconstellation(tLEO,dataLEO,GNSSFileType="YUMA");

    reg_tLEO(i) = tLEO;

    tLEO.Second=tLEO.Second+resolucion;

    Posicion_LEO(:,i) = satPos;

end

comet3(Posicion_LEO(1, :), Posicion_LEO(2, :), Posicion_LEO(3, :))

%%                              2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Hallar eventos de Radio Ocultación a partir de los dos criterios %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primer criterio
% eventos a una altitud de distancia mínima sobre la tierra entre -200 y
% 150 km

% Segundo criterio
% La ganancia de la antena se tienen en cuenta en el patrón de radiación


% Calculamos la linea de vista en cada punto de las trayectorias
% encontradas :
% Nota: Tener en cuenta que para una orbita completa de GPS aproximadamente
% 1 dia, tenemos varias orbitas LEO (Usat1) completas, a mayor velocidad.
% Es por ello que disminuir la resolución hace que tengamos más tiempo de
% procesamiento pero a su vez una mayor presición en la detección de
% eventos RO.

% Calculo de línea de vista y distancia minima al centro de la tierra
% Este procesamiento lo hacemos con SV-PRN(1) y la orbita LEO, luego
% lo repetimos para el resto de SV del almanaque.

clear reg_eventosRO
cont=1;
tam = size(Posicion_LEO);
hminLOS = zeros(1,length(tam(2)));

% Calculamos la altitud de la LOS que une punto por punto de las trayectorias.
% En realidad es la distancia minima de la recta al centro de la tierra

for index = 1:tam(2)

    l = Posicion_PRN1(:,index) - Posicion_LEO(:,index);
    p = -dot(l,Posicion_LEO(:,index))/norm(l)^2;
    minH = p*l+Posicion_LEO(:,index);
    minHlla = ecef2llaGeod(minH);
    hminLOS(index) = minHlla(3);

end

% Detectamos los cambios de signo en la altura

for i = 1: length(hminLOS)-1

    if hminLOS(i) > 0 && hminLOS(i+1) < 0

        % Si esto sucede tenemos un evento que tenia una altitud mayor y
        % paso a menor es un evento poniente

        % Guardamos los datos del evento y la posicion donde ocurrio para
        % despues ir a buscar el mometo de inicio y de final del evento

        % i es la posicion de GPS y LEO a la que se da el evento
        reg_eventosRO(cont,:) = [hminLOS(i) hminLOS(i+1) i] ;


        % VARIABLE QUE CONTIENE POS GPS Y LEO DEL EVENTO

        EVENTO_GPS(:,cont) = Posicion_PRN1(:,i);
        EVENTO_USAT(:,cont) = Posicion_LEO(:,i);
        cont = cont+1;
    elseif hminLOS(i) <0 && hminLOS(i+1) > 0

        % Si el evento sucede de negativo a positivo se registra un evento
        % de ocultación creciente
        % Esta variable guarda cual es la altitud (tiene que ser una positiva y una negativa)
        % ya que son los puntos donde cambia, detectamos en que
        reg_eventosRO(cont,:) = [hminLOS(i) hminLOS(i+1) i] ;

        % REG EVENTOS TIENE LA INFORMACION DE LA ALTITUD Y LA POSICION
        % DONDE SE DA EL EVENTO DE RO
        EVENTO_GPS(:,cont) = Posicion_PRN1(:,i);
        EVENTO_USAT(:,cont) = Posicion_LEO(:,i);

        cont = cont+1;

    end

end

% Situación gráfica
% tierra_color();
% hold on
%

% figure;
% Ground_track(Evento_RO(:,1:end-1)')

tierra_color()
hold on
for i=1:length(EVENTO_USAT)

    plot3([EVENTO_USAT(1,i) EVENTO_GPS(1,i)],[EVENTO_USAT(2,i) EVENTO_GPS(2,i)],[EVENTO_USAT(3,i) EVENTO_GPS(3,i)]);

    % pause(5)
end

 %%                                  3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Una vez que tenemos los momentos de positivo a negativo realizamos de
% nuevo una simulación pero con más resolución.

% cant event define la cantidad de eventos dentro del recorrido
% del satélite SVn y USAT que vamos a resolver
cant_event = length(reg_eventosRO);
% cant_event = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Posicion_USAT
clear Posicion_PRN

resolucion_n = 1; % [segundos]
ventana=300; % Tomamos una ventana de 400 segundos para la nueva propagación.
minHmin = [10e9;10e9;10e9]; % valiable para comparar la altitud, se setea este valor en alto para quitar la primera iteración
% Para cada uno de los anteriores eventos volvemos a simular en esta
% ventana.

for j=1:cant_event %length(reg_eventosRO)

    % Obtengo los tiempos de LEO y GPS de cada uno de los eventos
    minHmin = [10e9;10e9;10e9];
    reg_tGPS(reg_eventosRO(j,3)).Second = reg_tGPS(reg_eventosRO(j,3)).Second - 100; % Voy cargando el tiempo de inicio de los eventos
    tGPS_n = reg_tGPS(reg_eventosRO(j,3));

    reg_tLEO(reg_eventosRO(j,3)).Second = reg_tLEO(reg_eventosRO(j,3)).Second - 100; % Voy cargando el tiempo de inicio de los eventos
    tLEO_n = reg_tLEO(reg_eventosRO(j,3));


    % Empiezo 10 segundos antes para saber que el evento de RO está por
    % ocurrir

    clear reg_eventosRO_n
    clear hminLOS
    % Simulamos satélite y LEO

    for i=1:ventana

        % GPS

        [satPos] = gnssconstellation(tGPS_n,datos_PRN1,GNSSFileType="YUMA");

        tGPS_n.Second=tGPS_n.Second + resolucion_n;

        Posicion_PRN(:,i) = satPos;

        % LEO

        [satPos,satVel] = gnssconstellation(tLEO_n,dataLEO,GNSSFileType="YUMA");

        tLEO_n.Second=tLEO_n.Second + resolucion_n;

        Posicion_USAT(:,i) = satPos;
        Velocidad_USAT_f(:,i) = satVel;

    end
    
    %chequeo que el evento que analizamos acá sea correspondiente con el
    %que tenemos visto desde el código de arriba
    % tierra_color()
    % plot3([EVENTO_USAT(1,1) EVENTO_GPS(1,1)],[EVENTO_USAT(2,1) EVENTO_GPS(2,1)],[EVENTO_USAT(3,1) EVENTO_GPS(3,1)]);
    % comet3(Posicion_USAT(1, :), Posicion_USAT(2, :), Posicion_USAT(3, :))
    % break 
    tam = size(Posicion_USAT);

    cont=1;

    % Aplicamos el criterio de (+/-) 150 km a la LOS

    for index = 1:tam(2)

        l = Posicion_PRN(:,index) - Posicion_USAT(:,index);
        p = -dot(l,Posicion_USAT(:,index))/norm(l)^2;
        minH = p*l+Posicion_USAT(:,index);
        minHlla = ecef2llaGeod(minH);
        hminLOS(index) = minHlla(3);

        % Filtrado grueso a 10 seg para ver cuando hay cambio de signo
        % Tenemos el indice de posicion, la altura, hora, minuto, segundo.

        if (hminLOS(index) > -200e3 && hminLOS(index) < 150e3)
            reg_eventosRO_n(cont,:) = [index  hminLOS(index)];
            % Hacemos acá una detección de cual de los puntos de la
            % trayectoria tiene la altura minima para "definir" un evento
            a = ecef2llaGeod(minH);
            a = abs(a(3));
            b = ecef2llaGeod(minHmin);
            b = abs(b(3));
            if a < b

                minHmin = minH;

                % Definimos las posiciones de los satélites

                % Estos son los eventos de RO definidos

                EVENTO_USAT_fino(:,j)= Posicion_USAT(:,index);

                EVENTO_PRN_fino(:,j)= Posicion_PRN(:,index);

                EVENTO_PRN_vel(:,j) = Velocidad_USAT_f(:,index);

            end

            cont = cont + 1;

        end

    end


    % Nos armamos una variable que contenga la duración de cada evento y el
    % punto que "define" cada evento y la duración

    % [punto identificatorio del evento ; Duración del evento +-1seg ]


    Evento_RO(j,:) = [minHmin ; length(reg_eventosRO_n)];



    % Si queremos ver donde se mapea ese punto ejecutamos
    % Ground_track(Evento_RO(:,1:end-1)')

end







