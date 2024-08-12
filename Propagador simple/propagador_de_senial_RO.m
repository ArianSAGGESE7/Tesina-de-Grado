function [x,y,y0] = propagador_de_senial_RO(a)

% a es el numero de evento que queremos resolver (Event)
rLEO = load('RO_events_USAT.mat');
rLEO = rLEO.RO_events_USAT;
rGPS = load('RO_events_GPS.mat');
rGPS= rGPS.RO_events_GPS;
r = .1 ;
grad_step = 300; % Paso utilizado para el cálculo del gradiente del indice de refracción
step = 10000000; % Paso utilizado para el método de RKF (Runge Kutta Fehlberg)
exitos = 0; % Contador de pares GPS-LEO exitosos

q=1;
for index = a:a


    y0(:,1) = (rLEO(:,index)-rGPS(:,index))/norm(rLEO(:,index)-rGPS(:,index)); % Vector unitario en direccion del GPS al LEO
    y0(:,2) = llaGeod2ecef( ecef2llaGeod(rGPS(:,index)+y0(:,1)) + [0;0; -9.4979e-05 ]) - (rGPS(:,index)); % Vector unitario por encima de la elipsoide
    y0(:,2) = y0(:,2)/norm(y0(:,2));%

    [tau,x_prev,y_prev] =  RKF_aproximaciones(rGPS(:,index),rLEO(:,index),step,grad_step,y0(:,1),r); % Aplicamos una primera resolución con y0

    cnt = 0;

    while (1) % Bucle para los distintos disparos hasta converger con la precision deseada

        %   pause(5)

        [tau,x,y,n_reg,step_reg] = RKF_aproximaciones(rGPS(:,index),rLEO(:,index),step,grad_step,y0(:,cnt+2),r); % Resuelvo a partir de direcciones nuevas

        d = x(:,end) - rLEO(:,index); % "d" define la distancia entre el ultimo punto calculado y el punto donde yo queria llegar

        d_d(cnt+1) = norm(d); % Su norma debe establecer la toleracia "r" dada por nosotros

        if (norm (d) < r)
            exitos = exitos +1;
            ultima_posicion(:,exitos) = x(:,end); % esta para guardar la última posición del rayo y estudiar los angulos con el patrón de radiación de la antenna
            exitosos(exitos) = index;
            break;
        end

        % Corrección de la dirección a partir del método de la secante

        y0(:,cnt+3) = y0(:,cnt+2) + (rLEO(:,index) - x(:,end)) .* ((y0(:,cnt+2)-y0(:,cnt+1))./(x(:,end) - x_prev(:,end)));

        y0(:,cnt+3) = y0(:,cnt+3) / norm(y0(:,cnt+3));

        % Se reemplazan los valores viejos para estimar las direcciones
        % futuras
        x_prev = x;
        y_prev = y;

        % Limite de iteraciones: Si un par GPS-LEO tarda mas de 500
        % iteraciones en resolverse este se dará como fallido

        cnt = cnt + 1;
        if (cnt == 500)
            break;
        end

    end
    disparos_por_iteracion(index) = cnt;

    taus(index) = tau;
end

end