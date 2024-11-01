% Lazo completo -----------------------------------------------------------

clear all;clc;close

%--------------------------------------------------------------------------
%%                      Generación de señal sintética 
%--------------------------------------------------------------------------
NUMERO_DE_SATELITE = 1; % Adquisición del satélite 
cx = cacode (NUMERO_DE_SATELITE); % Adquiero un período de 1023 chips para el satélite elegido
fL1 = 1575.42e6; % Frecuencia nominal GPS
fOL = 1575e6; % Frecuencia de oscilador local (no es necesariamente la misma)
fFI = fL1-fOL; % Frecuencia intermedia
Tchip = 1/(1023e3); % Tiempo de chip nominal 
c = 3e8; % Velocidad de la luz
lambda = c/fL1; % Longitud de onda nominal 
fs = 5.3e6; % Frecuencia de muesteo predefinida
Ts = 1/fs; % Tiempo de muestreo 
fdata = 50; % Tasa de datos 
Tdata = 1/fdata; %Periodo de bit de datos (20ms)
CODE_LOOP          = 1; %Habilita lazo de código 
FREQ_LOOP          = 1; %Habilita lazo de portadora 
FC_ASIST           = 1; %Habilita asistencia e/ lazos 
%--------------------------------------------------------------------------
%                           Tiempo de simulación
TD = 2; % Duración de datos (seg)
n = 0:TD/Ts-1; % Indice de largo simulación
%--------------------------------------------------------------------------
%              Corrimiento por Doppler (consecuencia de retardos)
doppler = 100;
PEND = -doppler*lambda*Ts; % Como cambia el Doppler muestra a muestra
% Suponemos que para el tiempo de señal que queremos adquirir el receptor
% se mantiene cuasi estático con respecto al movimiento del satélite GPS
x = 20500e3 + (1:length(n))*PEND; % Rango [en metros] con inicialización en 20 km
% El fenómeno Doppler se traduce como un retardo temporal en las señales de
% banda base y portadora, por lo tanto existe un retardo asociado
taut  =x/c; % Tiempo asociado al pseudorango
% Cada señal se vé afectada por un pseudorango
cs = cx(mod(floor((n*Ts-taut)/Tchip),length(cx))+1);
ndata=0:TD/Tdata-1; %Indice para datos
data=sign(rand(1,length(ndata))-.5); % Datos generados de manera aleatoria
cdata=data(mod(floor((n*Ts-taut)/Tdata),length(data))+1);% Datos desplazados
% Con estás lineas de código literalmente se establece un retardo en el
% muestreo de los datos y los chips, ya que existe un Doppler que genera un
% efecto de shift en las muestras, aquí realizamos ese corrimiento
%--------------------------------------------------------------------------
%               Generación de la señal en envolvente compleja
s1 = cdata.*cs.*exp(1j*(2*pi*(fL1)*(Ts*n-taut))); % Señal modulada a fL1
s2 = s1.*exp(-1j*(2*pi*fOL*Ts*n)); % Demodulación a frecuencia intermedia
%--------------------------------------------------------------------------
%               Generación de ruido para mejorar el modelo 
CN0db = 50; % Relacion señal a ruido en DB
CN0 = 10^(.1*CN0db);
N = 1/(Ts*CN0);
% El ruido se genera a partir de una distribución complex normal
wI=randn(1,length(n));
wQ=randn(1,length(n));
nI=sqrt(N/2).*wI; %Ruido en fase
nQ=sqrt(N/2).*wQ; %Ruido en quadratura
ruido=nI+1i*nQ; %Ruido "Recibido"
%--------------------------------------------------------------------------
%                       Modelo de señal + ruido (sin considerar multicamino)
z = s2+ruido; 
Densidad_espectral(z,fs); 

%%                        Adquisición ( NO coherente)
%--------------------------------------------------------------------------
% Se define una cantidad de muestras a procesar en la adquisición 
Ti = 1e-3; % Tiempo de integración coherente
M = floor(Ti/Ts); % Cantidad de muestras a procesar
F = 10e3; % Frecuencias a recorrer
df = 0.1*(1/Ti); % Paso en frecuencias 
f_central = fFI; % Aquí está el espectro a la hora de adquirir
K = F/df ; % Cantidad de pasos en frecuencia 
feini  = f_central -F/2; % Frecuencia de inicio de la busqueda 
n = 1:M; % Ventana de integración
c=cx(mod(floor(n*Ts/Tchip),1023)+1); % Replica de código en banda base
cf=fft(c); %fft de código en banda base

INT_NC = 3; % Con esto indíco la cántidad de ventanas para la promediación
% se debe asegurar que existan esa cantidad de muestras sino se deberá
% correr desde más arriba una simulación más larga
% Se repite que acá no se consideran cambios de bit 
if(INT_NC*M<length(z))
    disp("Se corrió correctamete con la INT NC seleccionada");
end
rzt = 0;
for i=0:INT_NC-1
    fe = feini;
    for k = 1:K
        zp = z(i*M+1:(i+1)*M).*exp(-1j*(2*pi*fe*n*Ts)); % Quiero correlacionar con una replica
        % de código local y las señal que llegó
        zpf = fft(zp);
        rzs(k,:) = ifft(conj(cf).*zpf); % Función de intercorrelación
        fe = fe + df;
    end
    rzt = rzt + abs(rzs).^2;
end

rzs = rzt;
%% ------------------------------Graficos------------------------------------
figure
surf((1:M)*Ts/Tchip, -(feini:df:fe-1)+fFI,abs(rzs));
shading interp;
ylabel('$f-f_0$','Interpreter','latex');
xlabel('$\tau [chips]$','Interpreter','latex');
zlabel('$|r_{s\widetilde{s}}|$','Interpreter','latex');
title('Plano Retardo Doppler','Fontsize',14,'Interpreter','Latex');
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%%                             Busqueda del máximo
%--------------------------------------------------------------------------
% index contiene los indices de donde tengo que ir a buscar mi retardo de
% código y mi frecuencia 

[max_value, linear_index] = max(rzs(:));
[row, col] = ind2sub(size(rzs), linear_index);

a = (1:M)*Ts/Tchip;
b=(feini:df:fe)-fFI;

retardo_find =a(col); % Porcion en chips de retardo
frecuencia_find = b(row) +fFI; % Frecuencia de pico máximo (no centrada en fFI)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%                          Refinamiento en frecuencia
%--------------------------------------------------------------------------

% Se define una cantidad de muestras a procesar en la adquisición
Ti = 1e-3; % Tiempo de integración coherente
M = floor(Ti/Ts); % Cantidad de muestras a procesar
nt=1:M;
F = 200; % Frecuencias a recorrer
df_fino = 0.01*(1/Ti); % Cantidad de pasos en frecuencia
K = F/df_fino; % Cantidad de paso en frecuencia
feini = frecuencia_find-F/2;
dt_fino=.5; % Paso en retardo (en chips)
nc = 11; % Número de chips a recorrer
tauxini = retardo_find -(nc-1)/2;
zp = z(1:M); % Intervalo de datos a utilizar

fe = feini; % Frecuencia de inicio para la busqueda en el plano refinado
PC = 0; % Comparador para el pico

for k=1:M % Bucle para las frecuencias
    s = exp(1j*(2*pi*fe*nt*Ts));
    taux = tauxini;
    for i= 1:nc/dt_fino % Bucle para los retardos
        c=cx(mod(floor((nt*Ts-taux*Tchip)/Tchip),1023)+1); % Réplica de código
        rzs_fino(k,i) = sum(conj(c.*s).*zp); % Correlación discreta
        % (en vez de hacerlo con la FFT como el caso anterior)
        if (abs(rzs_fino(k,i))>PC)
            PC =abs(rzs_fino(k,i));
            fxe=fe;
            tauxe = taux;
        end
        taux = taux + dt_fino;
    end
    fe = fe +df_fino;
end

% Detectamos si el pico es lo suficientemente alto como para que sea
% declarado

UMB_co=4*sqrt(M*N); %Umbral de correlación


if(PC>UMB_co) %Si el candidato encontrado supera el umbral
    ADQ_OK=1;
    fx_fino=fxe; %Declaro frecuencia
    taux_fino=tauxe; %Declaro retardo (en chips)
end

%% Sin refinamiento: Pruebas
fx_fino=frecuencia_find ; %Declaro frecuencia
taux_fino=retardo_find;
%%
%------------------------------Graficos------------------------------------
figure
surf((1:M)*Ts/Tchip, -(feini:df:fe-1)+fFI,abs(rzs_fino));
shading interp;
ylabel('$f-f_0$','Interpreter','latex');
xlabel('$\tau [chips]$','Interpreter','latex');
zlabel('$|Rzs|$','Interpreter','latex');
title('Plano Retardo Doppler','Fontsize',14,'Interpreter','Latex');

%--------------------------------------------------------------------------
%%                      BUSQUEDA DE BIT DE DATOS
%--------------------------------------------------------------------------

T = 10e-3; % Tiempo de integración
M = floor(T/Ts); % Cantidad de muestras a procesar
 ntt=0:T/Ts-1; %Indice para ventana de integraci´ on de T segundos
 MSS=(TD/1e-3)-10; %Cantidad de correlaciones a hacer
 CDM=1e-3/Ts; %Cantidad de datos por milisegundo
% a una tasa de muestreo de fs, entran CMD datos.
% Replicas de portadora y de código
c=cx(mod(floor((ntt*Ts-taux_fino*Tchip)/Tchip),length(cx))+1);
s=exp(1j*(2*pi*fx_fino*ntt*Ts)); % no cambia la frecuencia por lo que 
% se actualiza una única vez.

 Tch_=1/(1/Tchip-1/(2*T)); %Tiempo de chip adulterado
 fdoppler=fx_fino-fFI;
 % Ayuda del lazo de portadora al de código, si no lo tenemos perfemos
 % potencia
 delta_T_frec=1/(fdoppler+fL1)-(1/fL1);
 delta_T_code=1023*1540*delta_T_frec;
 
for k=0:MSS-1
 zvd=z(CDM*k+1:CDM*(k+10)); %Selecci´ on de slot
 P_VD(k+1)=abs(sum(conj(c.*s).*zvd)); %Correlaci´ on
 c=cx(mod(floor((ntt*Ts-taux_fino*Tchip-(k)*delta_T_code)/Tchip),length(cx))+1); %R´ eplica de c´ odigo
 end

PROMPE=0.5*sum(P_VD)/length(P_VD);

plot(1:length(P_VD),P_VD);
% PRIM =0;
% for k=1:length(P_VD)
%     ANTERIOR=P_VD(k);
%     CENTRAL=P_VD(k+1);
%     POSTERIOR=P_VD(k+2);
%     if ((CENTRAL<POSTERIOR)&&(CENTRAL>ANTERIOR)&&(CENTRAL<PROMPE))
%                     MOM_TRAN=(k+4)*1e-3; %Instante en que ocurre la primera transición de bit de datos, tomando como referencia el origen de la ristra generada
%                     break;
%     end
% end
%--------------------------------------------------------------------------
%%                             Tracking
%--------------------------------------------------------------------------
MOM_TRAN = 89e-3;
fx=fx_fino; % Frecuencia (Estimación inicial)
taux=taux_fino+(fx-fFI)/1540*MOM_TRAN; % Retardo (Estimación inicial)
Ti = 10e-3; % Tomamos el tiempo de muestreo de 10 ms
M = floor(Ti/Ts) ; % Lo que nos da una cantidad de muestras a procesar de 
nt = 0:M-1; % Vecto de puntos que indican cuantas muestras son
MS = floor((TD-MOM_TRAN)/Ti); % Cantidad de slots de seguimiento
% Cuantos bloques de Ti segundos voy a realizar depende de la simulación

Tch_ = Tchip;

% Lazo de código 
BW_code = 1; % Ancho de banda del filtro en [Hz]
K0 = 4*BW_code*Ti; % Constante de lazo de código

% Lazo de portadora
BW_freq = 5; % BW de filtro de lazo de portadora
a2=2.4*BW_freq/0.7845*Ti; %Ctes de lazo de portadora
a1=1.1*(BW_freq/0.7845*Ti)^2;
a0=(BW_freq/0.7845*Ti)^3;

% Variables de estado del lazo de portadora 
out_p1=0;
out_p0=(fx-fFI)*2*pi*Ti;

c=cx(mod(floor((nt*Ts-taux*Tch_)/Tch_),length(cx))+1); %Réplica de código inicial
s=exp(1j*(2*pi*fx*nt*Ts)); %Réplica de portadora inicial


%Inicio de trackeo
for k=0:(MS-1)
    zseg=z(M*k+MOM_TRAN/Ts:M*(k+1)+MOM_TRAN/Ts-1); %Intervalo
    %Lazo de portadora:
    P=sum(conj(c.*s).*zseg); %Correlación prompt
    IPS=real(P);
    QPS=imag(P);
    D=atan(QPS/IPS); %Discriminador de fase
    out_p0=out_p0+a2*D; %variable de estado (Frecuencia estimada)
    out_p1=out_p1+a1*D+out_p0; %variable de estado (Fase estimada)
    fi_e=a0*D+out_p1; %Salida del lazo: fase estimada
    fdop=out_p0/(2*pi*Ti); %Frecuencia estimada

    %Lazo de código:
    sE=s.*cx(mod(floor((nt*Ts-taux*Tchip+.5*Tchip)/Tchip),1023)+1); %Réplica early
    E=sum(conj(sE).*zseg); %Correlación early
    sL=s.*cx(mod(floor((nt*Ts-taux*Tchip-0.5*Tchip)/Tchip),1023)+1); %Réplica late
    L=sum(conj(sL).*zseg); %Correlación late
    delta_tau=.5*(abs(E)-abs(L))/(abs(E)+abs(L)); %Discriminador
    f_tau=(fdop/1540); %Frecuencia del código alterada por doppler [chips/seg]

    taux=taux-K0*delta_tau+f_tau*Ti;  %Variable de estado del lazo de código y asistencia (realimentación)
    %del lazo de portadora
    s=exp(1j*(2*pi*(fFI+fdop)*Ts*nt+fi_e)); %Ajuste de réplica de portadora
    c=cx(mod(floor((nt*Ts-taux*Tch_)/Tch_),1023)+1); %Ajuste de réplica de código

    x_seg(k+1,1)=taux*Tchip*3e8; %Estimación de pseudorango
    fDD(k+1)=-fdop; %Estimación de doppler
    pE(k+1)=E; %Potencia de correlación early
    pL(k+1)=L; %Potencia de correlación late
    pP(k+1)=P; %Potencia de correlación prompt
    taux_vec(k+1)= taux;
end



close all
hold on


%Potencias de correlación
figure(1)
plot((0:length(pP)-1)*MS*Ti,abs(pP),'linewidth',1); grid on; hold on;
plot((0:length(pE)-1)*MS*Ti,abs(pE),'linewidth',1); hold on
plot((0:length(pL)-1)*MS*Ti,abs(pL),'linewidth',1); 

title(['Potencias P-E-L - SV',num2str(NUMERO_DE_SATELITE)],'Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
ylabel('Potencia [-]','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
legend('P','E','L')

%Estimación de Doppler
figure(2)
plot((0:length(pP)-1)*MS*Ti,fDD,'linewidth',1); grid on; hold on;
title('Doppler','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 

xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
legend('Doppler')

% Estimación retardo
figure(3)
plot((0:length(pP)-1)*MS*Ti,taux_vec,'linewidth',1); grid on; hold on;
title('Retardo','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 

xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
legend('Retardo')

%I y Q de potencia prompt
figure(4)
plot(1:10:10*length(pP),abs(real(pP)),'linewidth',1); grid on; hold on
plot(1:10:10*length(pP),imag(pP),'linewidth',1);
title(['Potencias I-Q Prompt - SV',num2str(NUMERO_DE_SATELITE)],'Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
ylabel('Potencia [-]','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
legend('En cuadratura ','En fase')
