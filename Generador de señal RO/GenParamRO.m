%========================================================================
%               Simulación de orbitas para extracto de parametros señales
%               de RO
%========================================================================
close all; clc ;clear variables

G = 6.67430e-11; % Constante de gravitación universal (m^3/kg/s^2)
M = 5.972e24;    % Masa de la Tierra (kg)
R = 6366000;     % Radio de la Tierra (m)


% Simulación de orbitas

data = yumaread("almanac.yuma.week0100.061440.txt");
SVn=3; % Satélite a simular
datos_PRN1 = data(SVn,:);
t = data.Time(1); % Tiempo asignado como inicial (estructura)
Periodo=2*pi*sqrt((datos_PRN1.SQRTA^2)^3/(M*G)); % Periodo de la órbita del satélite
resolucion = 1 ; % Cada cuanto tenemos una muestra
vueltas = 0.5;

cargar_datos =1;

if cargar_datos == 0

    for i=1:floor(Periodo/resolucion*vueltas)+resolucion
    
        [satPos,satVel] = gnssconstellation(t,datos_PRN1,GNSSFileType="YUMA");
    
        t.Second=t.Second+resolucion;
        tGPS(i) = t;
        PosSAT(:,i) = satPos';
        VelSAT(:,i) = satVel';
    end 
    % tierra_color()
    % comet3(PosSAT(1, :), PosSAT(2, :), PosSAT(3, :))
    
    % Simulación del satélite LEO
    dataLEO = yumaread("usat.yuma.txt");
    tLEO = dataLEO.Time(1);
    for i=1:(Periodo/resolucion)*vueltas+resolucion  
        [satPos,satVel] = gnssconstellation(tLEO,dataLEO,GNSSFileType="YUMA");    
        % tLEO(i) = tLEO;    
        tLEO.Second=tLEO.Second+resolucion;
        PosLEO(:,i) = satPos;
        VelLEO(:,i) = satVel;
    end
% comet3(PosLEO(1, :),PosLEO(2, :),PosLEO(3, :))
    else    
        datos = load("datosOrbit.mat");  
        VelLEO = datos.VelLEO;
        VelSAT = datos.VelSAT;
        PosLEO = datos.PosLEO;
        PosSAT = datos.PosSAT;
end

%% Calculo del Doppler 
VelLEO(:,1) = [VelLEO(1,2)/2 VelLEO(2,2) VelLEO(3,2)] ;
fL1 =1575.42e6;
c = 3e8;

rLOS = (PosLEO-PosSAT)./(vecnorm(PosLEO-PosSAT)); % Vector Linea de Vista
vRelLos = dot(VelLEO-VelSAT,rLOS); % Vector de proyección de la velocidad relativa 
fD = fL1/c * (vRelLos); % Doppler


%% Armado de la señal de radio Ocultación 

% Los pasos son los siguientes:

% Detectamos un evento con el criterio pre-establecido

% De ese vento medimos la altitud de la linea de vista y la grabamos 
% Esa linea de vista va a ser una altitud por una determinada cantidad de
% tiempo.


tam = size(PosLEO);
hminLOS = zeros(1,tam(2));
cont=1;

% Cálculo de la altitud de la LOS de cada satélite punto a punto

for index = 1:tam(2)

    l = PosSAT(:,index) - PosLEO(:,index);
    p = -dot(l,PosLEO(:,index))/norm(l)^2;
    minH = p*l+PosLEO(:,index);
    minHlla = ecef2llaGeod(minH);
    hminLOS(index) = minHlla(3);

end

% Determinación de eventos de RO (crecientes y decrecientes)

for i = 1: length(hminLOS)-1

    if hminLOS(i) > 0 && hminLOS(i+1) < 0

        PosGPS_Evento(:,cont) = PosSAT(:,i);
        PosLEO_Evento(:,cont) = PosLEO(:,i);
        cont = cont+1;
    elseif hminLOS(i) <0 && hminLOS(i+1) > 0

        PosGPS_Evento(:,cont) = PosSAT(:,i);
        PosLEO_Evento(:,cont) = PosLEO(:,i);
        cont = cont+1;

    end

end

cnt_row = 1;
cnt_col = 1;
flag =0;
for ii = 1:length(hminLOS)

    if (hminLOS(ii)>-150e3 && hminLOS(ii)<50e3)
    
            flag = 1;
            riestra_hminLOS(cnt_row,cnt_col) = hminLOS(ii);
            riestra_Doppler(cnt_row,cnt_col) = fD(ii);
            cnt_row = cnt_row +1;

    elseif (flag ==1)

            cnt_col = cnt_col +1;
            cnt_row =1;
            flag = 0;

    end

    

end


%% Interpolación LOS (Setea frecuencia de muestreo, tiempo de simulación)

% La idea es poder generar una funcion de LOS que me permita agarrar
% distintos eventos que estan en la riestra y poder interpolarlos para
% muestrearlos nuevamente. (tener en cuenta que las muestras estan 1 por lo que indica la variable 'resolución')

evento = 5; % Selecciono cuales de los eventos en donde la LOS cambia de signo (o sea hubo un evento).

MuestrasLOSEvento = riestra_hminLOS(:,evento);
MuestrasLOSEvento = MuestrasLOSEvento(MuestrasLOSEvento ~= 0); % Sacamos los ceros porque no todos los eventos duran lo mismo (y se rellena por default con cero)

MuestrasDopplerEvento = riestra_Doppler(:,evento); % Tomo de ese evento que registro la alitud del punto minimo de LOS su respectivo Doppler
MuestrasDopplerEvento = MuestrasDopplerEvento(MuestrasDopplerEvento~= 0);

tEvento = (0:length(MuestrasLOSEvento)-1)*resolucion; % Esto es lo que dura el evento que analizamos

fs = 2.3e6; % Tasa de muestreo ====================================
Ts = 1/fs;
Tsim = 100;
tInterp= (0:Ts:Tsim); % Tiempo a evaluar (con la resolución que se requiera, esto lo fijamos nosotros)
N = length(tInterp); % Cantidad de muestras a procesar 


%% Generación de señal en banda base solo Doppler

MuestrasLOSInterp = interp1(tEvento,MuestrasLOSEvento, tInterp, 'spline'); % LOS para esa resolución
MuestrasDopplerInterp = interp1(tEvento,MuestrasDopplerEvento,tInterp,'spline'); % Doppler para esa resolución 


% Se tiene que interpolar tanto la amplitud como la fase de forma tal de
% que sucedan para lo que marca tInterp

% Interpolación de la amplitud en función de la hLOS
amplitude_LEO = load('datosGenRO.mat', 'amplitude_LEO');
phase_LEO = load('datosGenRO.mat', 'phase_LEO');

h_LOS_datos = linspace(-524.288e3/2, 524.288e3/2, 2^10)'; % Datos que pertenecen a las muestras reales (simulados)

amplitude_LEO = amplitude_LEO.amplitude_LEO'; % Datos que pertenecen a las muestras reales (simulados)
phase_LEO = phase_LEO.phase_LEO';


Amplitude_AJ = interp1(h_LOS_datos, amplitude_LEO, MuestrasLOSInterp, 'spline'); % Interpolamos amplitud
Phase_INT = interp1(h_LOS_datos, phase_LEO, MuestrasLOSInterp, 'spline'); % Interpolamos la fase

% Cada evento tiene una duración máxima intrínseca de la duración real,
% nosotros interpolamos los parametros los cuales tienen como rango de
% validez hasta la duración real.

%==========================================================================
%                       Recreamos nuestra señal de RO
%==========================================================================
% El tiempo queda definido por Tsim, para ese tiempo extraemos dopler y una
% dada hLOS

fL1 = 1545.75e6;
fL1 = 0;  % Lo dejamos así para seguir unicamente las variaciones Doppler
ampRO = Amplitude_AJ;
phaseRO = Phase_INT;
fD_GEOM = MuestrasDopplerInterp;

sRO = ampRO.*exp(1j*(2*pi*cumtrapz(tInterp,(fL1 + fD_GEOM))+ phaseRO)); % con fase del evento
% sRO = ampRO.*exp(1j*(2*pi*cumtrapz(tInterp,(fL1 + fD_GEOM))); % sin fase del evento 
% sRO = exp(1j*(2*pi*cumtrapz(tInterp,(fL1 + fD_GEOM)))); % sin fase del evento sin modificación de la amplitud 

% Se le debe aplicar ruido a la señal para simular el entorno, si queremos
% una CN0 de 45 db al final del evento tomamos 
A = 0.8 ; % Amplitud de la portadora para una altitud de 0 km de LOS 
CN0_db = 45;
CN0 = 10^(0.1*CN0_db);
N0 = A^2/2/(CN0);
N0_var = fs*N0;
 
% El ruido se genera a partir de una distribución complex normal
wI=randn(1,N);
wQ=randn(1,N);
nI=sqrt(N0_var/2).*wI; %Ruido en fase
nQ=sqrt(N0_var/2).*wQ; %Ruido en quadratura
ruido=nI+1i*nQ; %Ruido "Recibido"


sRO1 = sRO + ruido; % señal de salida para guardar 

% Guardar simulación
etiqueta = sprintf('Frecuencia de muestreo: %d Hz\nNúmero de muestras: %d \n CN0 [dB] = %d', fs, N,CN0_db);
datosSRO.muestras = sRO1;
datosSRO.etiqueta = etiqueta; 
datosSRO.doppler = fD_GEOM; 
datosSRO.amplitud = ampRO;
datosSRO.fase = phaseRO;

if MuestrasLOSInterp(1) < 0
    datosSRO.evento = 1; % Creciente (se invierte)

else 
    datosSRO.evento = 0; % Decreciente (no se invierte)
end
save('Señal_RO', 'datosSRO', '-v7.3','-nocompression');
%% Generación señal GNSS (chips + datos)

% Interpolación de los datos preexistentes

MuestrasLOSInterp = interp1(tEvento,MuestrasLOSEvento, tInterp, 'spline'); % LOS para esa resolución
MuestrasDopplerInterp = interp1(tEvento,MuestrasDopplerEvento,tInterp,'spline'); % Doppler para esa resolución 


% Se tiene que interpolar tanto la amplitud como la fase de forma tal de
% que sucedan para lo que marca tInterp

% Interpolación de la amplitud en función de la hLOS
amplitude_LEO = load('datosGenRO.mat', 'amplitude_LEO');
phase_LEO = load('datosGenRO.mat', 'phase_LEO');

h_LOS_datos = linspace(-524.288e3/2, 524.288e3/2, 2^10)'; % Datos que pertenecen a las muestras reales (simulados)

amplitude_LEO = amplitude_LEO.amplitude_LEO'; % Datos que pertenecen a las muestras reales (simulados)
phase_LEO = phase_LEO.phase_LEO';


Amplitude_AJ = interp1(h_LOS_datos, amplitude_LEO, MuestrasLOSInterp, 'spline'); % Interpolamos amplitud
Phase_INT = interp1(h_LOS_datos, phase_LEO, MuestrasLOSInterp, 'spline'); % Interpolamos la fase

fL1 = 1545.75e6;
fL1 = 0;  % Lo dejamos así para seguir unicamente las variaciones Doppler
ampRO = Amplitude_AJ;
phaseRO = Phase_INT;
fD_GEOM = MuestrasDopplerInterp;


% Datos de señal GNSS-GPS

NUMERO_DE_SATELITE = 1; % Adquisición del satélite 
cx = cacode (NUMERO_DE_SATELITE); % Adquiero un período de 1023 chips para el satélite elegido
fL1 = 1575.42e6; % Frecuencia nominal GPS
fOL = 1575e6; % Frecuencia de oscilador local (no es necesariamente la misma)
fFI = fL1-fOL; % Frecuencia intermedia
Tchip = 1/(1023e3); % Tiempo de chip nominal 
C = 3e8; % Velocidad de la luz
fdata = 50; % Tasa de datos 
Tdata = 1/fdata; %Periodo de bit de datos (20ms)

taut = 0; % Este es el retardo que tenemos que modelar 

cs = cx(mod(floor((tInterp-taut)/Tchip),length(cx))+1); % Código
ndata=0:TD/Tdata-1; % Indice de datos
data=sign(rand(1,length(ndata))-.5); % Datos generados de manera aleatoria
cdata=data(mod(floor((tInterp-taut)/Tdata),length(data))+1);% Datos desplazados

sRO = ampRO.*exp(1j*(2*pi*cumtrapz(tInterp,(fD_GEOM))+ phaseRO)); % Señal en banda base sin retardo

A = 0.8 ; % Amplitud de la portadora para una altitud de 0 km de LOS 
CN0_db = 45;
CN0 = 10^(0.1*CN0_db);
N0 = A^2/2/(CN0);
N0_var = fs*N0;
 
% El ruido se genera a partir de una distribución complex normal
wI=randn(1,N);
wQ=randn(1,N);
nI=sqrt(N0_var/2).*wI; %Ruido en fase
nQ=sqrt(N0_var/2).*wQ; %Ruido en quadratura
ruido=nI+1i*nQ; %Ruido "Recibido"

sROGNSS = sRO + ruido;
saveVector(sROGNSS, 1e6, SigRoGnss);
%% Graficos de interpolación y ajuste
close all
figure;
subplot(4,1,1)
hold on
plot(tEvento,MuestrasLOSEvento,'r','LineWidth',1.5);
plot(tInterp,MuestrasLOSInterp,'b','LineWidth',1);
title('Interpolación LOS')
subplot(4,1,2)
hold on
plot(tEvento,MuestrasDopplerEvento,'r','LineWidth',1.5);
plot(tInterp,MuestrasDopplerInterp,'b','LineWidth',1);
title('Doppler')
subplot(4,1,3)
hold on 
plot(h_LOS_datos,amplitude_LEO,'r','LineWidth',1.5);
plot(MuestrasLOSInterp,Amplitude_AJ,'b','LineWidth',1);
title('Ajuste de amplitud en función de LOS')
subplot(4,1,4)
hold on 
plot(h_LOS_datos,phase_LEO,'r','LineWidth',1.5);
plot(MuestrasLOSInterp,Phase_INT,'b','LineWidth',1);
title('Interpolación de fase en función en función de LOS')

figure;
plot(tInterp,real(sRO),'b','LineWidth',1)
title('Parte Real de la señal de RO')
% Tener en cuenta que el eje x en los últimos dos graficos corresponde a
% LOS y no a tiempo ya que de ahi sacamos las simulaciones anteriores
%% Graficos de Doppler de simulación de órbita
figure;
subplot(3,1,1)
plot(tInterp,fD,'LineWidth',1); 
ylabel('$f_{Doppler}$','Interpreter','latex')
subplot(3,1,2)
plot((0:length(hminLOS)-1)*resolucion,hminLOS,'LineWidth',1)
ylabel('$h_{LOS}$','Interpreter','latex')
subplot(3,1,3)
plot(diff(hminLOS),'LineWidth',1)
ylabel('$dh_{LOS} / dt$','Interpreter','latex')
xlabel('tiempo [Segundos]','Interpreter','latex')