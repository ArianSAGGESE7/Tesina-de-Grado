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
SVn = 1; % Satélite a simular
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

d = vecnorm(PosLEO-PosSAT);

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
            riestra_distancia(cnt_row,cnt_col) = d(ii); 
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

evento = 8; % Selecciono cuales de los eventos en donde la LOS cambia de signo (o sea hubo un evento).

MuestrasLOSEvento = riestra_hminLOS(:,evento);
MuestrasLOSEvento = MuestrasLOSEvento(MuestrasLOSEvento ~= 0); % Sacamos los ceros porque no todos los eventos duran lo mismo (y se rellena por default con cero)


MuestrasDopplerEvento = riestra_Doppler(:,evento); % Tomo de ese evento que registro la alitud del punto minimo de LOS su respectivo Doppler
MuestrasDopplerEvento = MuestrasDopplerEvento(MuestrasDopplerEvento~= 0);


MuestrasDistancia = riestra_distancia(:,evento);
MuestrasDistancia = MuestrasDistancia(MuestrasDistancia~= 0);

tEvento = (0:length(MuestrasLOSEvento)-1)*resolucion; % Esto es lo que dura el evento que analizamos
fs = 105001; % track tono
% fs = 2092000; % Tasa de muestreo ====================================     CAMBIAR ESTO AL SCRIPT DONDE VAYAS A CORRER
Ts = 1/fs;
Tin = 0;
Tsim = 85;
tInterp= (Tin:Ts:Tsim); % Tiempo a evaluar (con la resolución que se requiera, esto lo fijamos nosotros)
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
datosSRO.fs = fs;
datosSRO.CN0 = CN0_db;
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
MuestrasDistanciaInterp = interp1(tEvento,MuestrasDistancia,tInterp,'spline'); % Distancia entre satelites para esa resolución 

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

ampRO = Amplitude_AJ;
phaseRO = Phase_INT;
fD_GEOM = MuestrasDopplerInterp;


% Datos de señal GNSS-GPS

NUMERO_DE_SATELITE = 1; % Adquisición del satélite 
cx = cacode (NUMERO_DE_SATELITE); % Adquiero un período de 1023 chips para el satélite elegido
fL1 = 1575.42e6; % Frecuencia nominal GPS
fOL = 1575.42e6; % Frecuencia de oscilador local (no es necesariamente la misma)
fFI = fL1-fOL; % Frecuencia intermedia
Tchip = 1/(1023e3); % Tiempo de chip nominal 
C = 3e8; % Velocidad de la luz
fdata = 50; % Tasa de datos 
Tdata = 1/fdata; %Periodo de bit de datos (20ms)

% Calculo del retardo
Vrel = C*fD_GEOM/fL1;
dis = 20000e3+ cumtrapz(tInterp,(Vrel)); % Modelo del retardo de la señal
taut = dis/C;

% Calculo retardo con otra pend Doppler
% doppler = 35560;
% PEND = -doppler.*C/fL1*Ts;
% x = 20000e3 + (1:length(tInterp)).*PEND;
% taut = x/C;

cs = cx(mod(floor((tInterp-taut)/Tchip),length(cx))+1); % Código
ndata=0:Tsim/Tdata-1; % Indice de datos
data=sign(rand(1,length(ndata))-.5); % Datos generados de manera aleatoria
cdata=data(mod(floor((tInterp-taut)/Tdata),length(data))+1);% Datos desplazados


sRO = cdata.*cs.*ampRO.*exp(1j*(2*pi*-taut*fL1+ phaseRO))/ampRO(1); % Señal en banda base sin retardo
% sRO = cdata.*cs.*exp(1j*(2*pi*-1*taut*fL1)); % Señal GNSS pura en banda base

% Recordar que la señal tiene que ser muestreada con fs maor a 2M por el
% código
CN0db = 45; % Relacion señal a ruido en DB
CN0 = 10^(.1*CN0db);
N = 1/(Ts*CN0);
% El ruido se genera a partir de una distribución complex normal
wI=randn(1,length(sRO));
wQ=randn(1,length(sRO));
nI=sqrt(N/2).*wI; %Ruido en fase
nQ=sqrt(N/2).*wQ; %Ruido en quadratura
ruido=nI+1i*nQ; %Ruido "Recibido"

sROGNSS = sRO + ruido;
% saveVector(sROGNSS, 50e6, 'SigRoGnss');

%%
% Guardar simulación
etiqueta = sprintf('Frecuencia de muestreo: %d Hz\nNúmero de muestras: %d \n CN0 [dB] = %d', fs, N,CN0db);
datosSRO.etiqueta = etiqueta;
datosSRO.muestras = sROGNSS; 
datosSRO.doppler = fD_GEOM; 
datosSRO.amplitud = ampRO;
datosSRO.fase = phaseRO;
datosSRO.fs = fs;
datosSRO.CN0 = CN0db;

if MuestrasLOSInterp(1) < 0
    datosSRO.evento = 1; % Creciente (se invierte)
else 
    datosSRO.evento = 0; % Decreciente (no se invierte)
end
save('Señal_RO', 'datosSRO', '-v7.3','-nocompression');
%% Graficos de interpolación y ajuste
close all
grayColor =[0.4940 0.1840 0.5560];
figure;
subplot(4,1,1)
hold on 
box on
grid on
plot(tEvento,MuestrasLOSEvento,'r','LineWidth',1.5);
plot(tInterp,MuestrasLOSInterp,'Color',grayColor,'LineWidth',2);
title('$h_{min_{LOS}}$','Interpreter','latex')
xlabel('tiempo [seg]','Interpreter','latex')
subplot(4,1,2)
hold on 
box on
grid on
plot(tEvento,MuestrasDopplerEvento,'r','LineWidth',2);
plot(tInterp,MuestrasDopplerInterp,'Color',grayColor,'LineWidth',2);
title('$f_{DG}$','Interpreter','latex')
xlabel('tiempo [seg]','Interpreter','latex')
subplot(4,1,3)
hold on 
box on
grid on
plot(h_LOS_datos,amplitude_LEO,'k','LineWidth',2);
plot(MuestrasLOSInterp,Amplitude_AJ,'Color',grayColor,'LineWidth',2);
ylabel('$\hat{A}_{RO}$','Interpreter','latex')
title('Interpolaci\''on de amplitud en funci\''on de LOS','Interpreter','latex')
xlabel('$h_{min_{LOS}}$','Interpreter','latex')
subplot(4,1,4)
hold on 
box on
grid on
plot(h_LOS_datos,phase_LEO,'k','LineWidth',2);
plot(MuestrasLOSInterp,Phase_INT,'Color',grayColor,'LineWidth',2);
ylabel('$\hat{\Phi}_{E}$','Interpreter','latex')
title('Interpolaci\'' on de fase en funci\'' on de LOS','Interpreter','latex')
xlabel('$h_{min_{LOS}}$','Interpreter','latex')

figure
hold on
box on
grid on
plot(tInterp,real(sRO),'k','LineWidth',2)
title('Parte Real de la $se\tilde{n}al$ de RO','Interpreter','latex')
xlabel('tiempo [seg]','Interpreter','latex')
ylabel('Amplitud','Interpreter','latex')
% Tener en cuenta que el eje x en los últimos dos graficos corresponde a
% LOS y no a tiempo ya que de ahi sacamos las simulaciones anteriores
%% Graficos de Doppler de simulación de órbita
figure;
subplot(3,1,1)
    hold on
    plot((0:length(hminLOS)-1)*resolucion,hminLOS,'k','LineWidth',2)
    ylabel('$h_{min_{LOS}}$','Interpreter','latex')
    box on
    grid on
    xlim([0 21538])
    
for i =1: length(hminLOS)
    
    if (hminLOS(i)>-150e3 && hminLOS(i)<50e3)
    plot((i-1)*resolucion,hminLOS(i),'-s','MarkerSize',5,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[0 0.4470 0.7410])
    end
end
title('Simulaci\'' on para generaci\'' on de $se\tilde{n}al$ GNSS-RO','Interpreter','latex')
xlabel('a)','Interpreter','latex')
subplot(3,1,2)

hold on
for i =1: length(hminLOS)
    
    if (hminLOS(i)>-150e3 && hminLOS(i)<50e3)
    plot((i-1)*resolucion,fD(i),'-s','MarkerSize',5,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[0 0.4470 0.7410])
    end
end


plot((0:length(fD)-1)*resolucion,fD,'k','LineWidth',2); 
ylabel('$f_{DG}$','Interpreter','latex')
box on
grid on
xlim([0 21538])
xlabel('b)','Interpreter','latex')
subplot(3,1,3)

plot(diff(hminLOS),'k','LineWidth',2)
ylabel('$dh_{LOS} / dt$','Interpreter','latex')
xlabel('c) tiempo [Segundos]','Interpreter','latex')
box on
grid on
xlim([0 21538])