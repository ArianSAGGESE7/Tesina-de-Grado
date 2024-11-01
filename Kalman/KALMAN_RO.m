%==========================================================================
% Implementación del Filtro de Kalman para mediciónes de RADIO OCULTACIÓN 
%==========================================================================
                        clc;clear all;close all;


%               ===========================================
%             Primera prueba con señal sintética normal de GPS
%               ===========================================
%% Carga de datos 

NUMERO_DE_SATELITE = 1; % Adquisición del satélite 
cx = cacode (NUMERO_DE_SATELITE); % Adquiero un período de 1023 chips para el satélite elegido
fL1 = 1575.42e6; % Frecuencia nominal GPS
fOL = 1575e6; % Frecuencia de oscilador local (no es necesariamente la misma)
fFI = fL1-fOL; % Frecuencia intermedia
Tchip = 1/(1023e3); % Tiempo de chip nominal 
C = 3e8; % Velocidad de la luz
lambda = C/fL1; % Longitud de onda nominal 
fs = 5.3e6; % Frecuencia de muesteo predefinida
Ts = 1/fs; % Tiempo de muestreo 
fdata = 50; % Tasa de datos 
Tdata = 1/fdata; %Periodo de bit de datos (20ms)


% Enablear gráficos
fig_adq = 1; % Graficos de la sección de adquisición

%% Generación de señal de GPS sintética 

TD = 6; % Duración de datos (seg)
n = 0:TD/Ts-1; % Indice de largo simulación
doppler = 100.3333;
PEND = -doppler*lambda*Ts; % Como cambia el Doppler muestra a muestra
% Suponemos que para el tiempo de señal que queremos adquirir el receptor
% se mantiene cuasi estático con respecto al movimiento del satélite GPS
x = 20500e3 + (1:length(n))*PEND; % Rango [en metros] con inicialización en 20 km
% El fenómeno Doppler se traduce como un retardo temporal en las señales de
% banda base y portadora, por lo tanto existe un retardo asociado
taut  =x/C; % Tiempo asociado al pseudorango
% Cada señal se vé afectada por un pseudorango
cs = cx(mod(floor((n*Ts-taut)/Tchip),length(cx))+1);
ndata=0:TD/Tdata-1; %Indice para datos
data=sign(rand(1,length(ndata))-.5); % Datos generados de manera aleatoria
cdata=data(mod(floor((n*Ts-taut)/Tdata),length(data))+1);% Datos desplazados
s1 = cdata.*cs.*exp(1j*(2*pi*(fL1)*(Ts*n-taut))); % Señal modulada a fL1
s2 = s1.*exp(-1j*(2*pi*fOL*Ts*n)); % Demodulación a frecuencia intermedia
%Generación de ruido para mejorar el modelo 
CN0db = 40; % Relacion señal a ruido en DB
CN0 = 10^(.1*CN0db);
N = 1/(Ts*CN0);
% El ruido se genera a partir de una distribución complex normal
wI=randn(1,length(n));
wQ=randn(1,length(n));
nI=sqrt(N/2).*wI; %Ruido en fase
nQ=sqrt(N/2).*wQ; %Ruido en quadratura
ruido=nI+1i*nQ; %Ruido "Recibido"
z = s2+ruido; % Modelo de señal + ruido 

%% Etapa de ADQUISICIÓN

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

rzs=rzt;
% Busqueda del máximo 
[max_value, linear_index] = max(rzs(:));
[row, col] = ind2sub(size(rzs), linear_index);

a = (1:M)*Ts/Tchip;
b=(feini:df:fe)-fFI;
if (b(row) == 0)

    row = row+1;
end
retardo_find =a(col); % Porcion en chips de retardo
frecuencia_find = b(row) +fFI; % Frecuencia de pico máximo (no centrada en fFI)

if fig_adq

    figure
    surf((1:M)*Ts/Tchip, -(feini:df:fe-1)+fFI,abs(rzs));
    shading interp;
    ylabel('$f-f_0$','Interpreter','latex');
    xlabel('$\tau [chips]$','Interpreter','latex');
    zlabel('$|r_{s\widetilde{s}}|$','Interpreter','latex');
    title('Plano Retardo Doppler','Fontsize',14,'Interpreter','Latex');

end
%% Identificamos el comienzo del bit de datos 

fx_fino=frecuencia_find ; %Declaro frecuencia
taux_fino=retardo_find;
T = 10e-3; % Tiempo de integración
M = floor(T/Ts); % Cantidad de muestras a procesar
ntt=0:T/Ts-1; %Indice para ventana de integración de T segundos
MSS=(TD/1e-3)-10; %Cantidad de correlaciones a hacer
CDM=1e-3/Ts; %Cantidad de datos por milisegundo a una tasa de muestreo de fs, entran CMD datos.
c=cx(mod(floor((ntt*Ts-taux_fino*Tchip)/Tchip),length(cx))+1);% Replicas de portadora y de código
s=exp(1j*(2*pi*fx_fino*ntt*Ts)); % no cambia la frecuencia por lo que se actualiza una única vez.
fdoppler=fx_fino-fFI;
delta_T_frec=1/(fdoppler+fL1)-(1/fL1);% Ayuda del lazo de portadora al de código, si no lo tenemos perfemos% potencia
delta_T_code=1023*1540*delta_T_frec;
 
for k=0:MSS-1
 zvd=z(CDM*k+1:CDM*(k+10)); %Selección de slot
 P_VD(k+1)=abs(sum(conj(c.*s).*zvd)); %Correlación
 c=cx(mod(floor((ntt*Ts-taux_fino*Tchip-(k)*delta_T_code)/Tchip),length(cx))+1); %Réplica de código
end

PROMPE=0.5*sum(P_VD)/length(P_VD);

plot(1:length(P_VD),P_VD); %con este gráfico podemos identificar a
% partir de cuando se puede empezar con los lazos y el filtrado para no
% incurrir en un cambio de signo por parte de los bit de datos
%% Etapa de inicialización del Filtro y Loop 

Ti = 10e-3; % Tomamos el tiempo de integración de 5 ms
T=Ti;
% Parámetros del FIltro
F= [1 T T^2/2; 0 1 T; 0 0 1]; % Matriz de transición de estados 
H = [1 T/2 T^2/6]; % Matriz de medición 


% H = [1 0 0]; % Matriz de medición 


FL =  fL1; % Frecuencia de portadora GPS
h0 = 2e-23; % Parámetros de ruido
h2 = 2.13e-24;
q0 = h0/2;
qw = 2*pi^2*h2;
qa = 0*10^(-10); % Se obtiene de una tabla que se muestra en libro de Montenbruk;
c = 3e8;

% Matriz de diseño para tracking

Q = qa*[T^5/20 T^4/8 T^3/6; T^4/8 T^3/3 T^2/2; T^3/6 T^2/2 T] + q0*[T^3/3 T^2/2 0; T^2/2 T 0; 0 0 0] + qw*[T 0 0; 0 0 0; 0 0 0];


% Matriz paper RO after tracking

% Q = 1e6*[(T*q0+T^3/3*qw + T^5/20*qa/c^2)  (T^2/2*qw + T^4/8*qa/c^2)  (T^3/6*qa/c^2);
%      (T^2/2*qw + T^4/8*qa/c^2)  (T*qw + T^3/3*qa/c^2) (T^2/2*qa/c^2)
%      (T^3/6*qa/c^2)  (T^2/2*qa/c^2) (T*qa/c^2)]*(2*pi*FL)^2; % Covarianza del proceso de ruido 


% Matriz de Fermin

% h_1=0;
% h_2=2.13e-24;
% Sj = 90e3;   % Intensidad del jerk. "Smoother-based GPS Signal Tracking in a SDR" pp. 22 da la refe de Sj=1300 rad^2/s^5
% Q = 5e-4*[pi^2*T^5/5 pi*T^4/4 pi*T^3/3;pi*T^4/4 T^3/3 T^2/2;pi*T^3/3 T^2/2 T]*Sj;
% % Allan variance covariance model (2-state)
% covxy = (2*pi*fL1)^2*[(h0*T + 2*h_1*T^2 + 2/3*pi^2*h_2*T^3) , (2*h_1*T + pi^2*h_2*T^2) ; ...
%                               (2*h_1*T + pi^2*h_2*T^2)              , (2*h_1 + 2*pi^2*h_2*T)];
% Q = Q + [covxy, [0 0].';[0 0 0]];


P0 = [pi^2/3 0 0; 0 (2*pi*100)^2/12 0 ; 0 0 5] ; % Matriz de covarianza del error (Hay que incializarla)

x0 = [0;2*pi*(fdoppler);0]; % Inicialización de los estados con lo que obtuvimos de la etapa de ADQ


% Parámetros del filtro de lazo
BW_code = 1; % Ancho de banda del filtro en [Hz]
K0 = 4*BW_code*Ti; % Constante de lazo de código

close all
MOM_TRAN = 28e-3; % Esta variable define el tiempo a a partir del cual puedo integrar
fx=fx_fino; % Frecuencia (Estimación inicial)
taux=taux_fino+(fx-fFI)/1540*MOM_TRAN; % Retardo (Estimación inicial)
M = floor(Ti/Ts) ; % Lo que nos da una cantidad de muestras a procesar de 
nt = 0:M-1; % Vecto de puntos que indican cuantas muestras son
MS = floor((TD-MOM_TRAN)/Ti); % Cantidad de slots de seguimiento

x = zeros(3,MS);
Kk = zeros(3,1,MS);    
Px_post = zeros(3,3,MS);
Px_prior = zeros(3,3,MS);
Px_post(:,:,1) = P0;
x(:,1)= x0;


% Construyo las réplicas iniciales de la señal compleja
c=cx(mod(floor((nt*Ts-taux*Tchip)/Tchip),length(cx))+1); %Réplica de código inicial
s=exp(1j*(2*pi*fx*nt*Ts)); %Réplica de portadora inicial

% Se interpreta en que se toman de a slots de M muestras de un total de MS muestras, entre cálculo y
% cálculo pasaron, cada instante en muestras equivale a un Ti

for k=0:(MS-2)
    zseg=z(M*k+MOM_TRAN/Ts:M*(k+1)+MOM_TRAN/Ts-1); % Intervalo de muestras 

    % Correladores Prompt, Early y Late
    P=sum(conj(c.*s).*zseg); %Correlación prompt
    sE=s.*cx(mod(floor((nt*Ts-taux*Tchip+.5*Tchip)/Tchip),1023)+1); %Réplica early
    E=sum(conj(sE).*zseg); %Correlación early
    sL=s.*cx(mod(floor((nt*Ts-taux*Tchip-0.5*Tchip)/Tchip),1023)+1); %Réplica late
    L=sum(conj(sL).*zseg); %Correlación late

    % Discriminadores
    Ip=real(P);
    Qp=imag(P);
    % D_fase = atan2(Qp,Ip); %Discriminador de fase
    D_fase = atan(Qp/Ip);
    D_tau=.5*(abs(E)-abs(L))/(abs(E)+abs(L)); %Discriminador de código
    
    % Componentes en fase y cuadratura para el cálculo de la CN0
    IQ = [real(zseg);imag(zseg)]; % primera fila datos en fase y segunda en cuadratura en Ti segundos de integración

    % Estimación del retardo
    taux=taux-K0*D_tau+x(2,k+2)*Ti/1540/2/pi; 

    % Estimación por Kalman
    est_x_prior = F*x(:,k+1); % Propagación
    Px_prior(:,:,k+2) = F*Px_post(:,:,i+1)*F.' + Q; % Estimación
    % Para el cálculo de R debemos estimar CN0
    % [CN0]=CN0_estimation(IQ(1,:),IQ(2,:),Ti);
    % CN0_vec(k+1) = CN0;
    % CN0=10^(0.1*40);
    R = 1/(2*Ti*CN0); % Podriamos acá agregar la parte que esta en el paper con SI
    K = Px_prior(:,:,k+2)*H.'/(H*Px_prior(:,:,k+2)*H.' + R);% Ganancia de Kalman

    est_x_posterior = est_x_prior + K*(D_fase);

    Px_post(:,:,k+2) = Px_prior(:,:,k+2) - K*H*Px_prior(:,:,k+2); % Cálculo de innovaciones
    
    % Vectores
    x(:,k+2) = est_x_posterior;
    Kk(:,:,k+2) = K; 
    Pp(:,k+1)=P;
    Ep(:,k+1)=E;
    Lp(:,k+1)=L;
    
    % Generamos réplica local para volver a entrar a los discriminadores 
    s=exp(1j*(2*pi*(fFI+x(2,k+2)/2/pi)*Ts*nt+x(1,k+2))); % Ajuste de réplica de portadora
    c=cx(mod(floor((nt*Ts-taux*Tchip)/Tchip),1023)+1); % Ajuste de réplica de código

    % Continuamos en el lazo
end

%% Gráficos
close all
hold on 
% figure(1)
% plot((0:length(x(1,:))-1)*Ti,x(1,:),'LineWidth',1) % Fase
figure(1)
plot((0:length(x(1,:))-1)*Ti,x(2,:)/2/pi,'LineWidth',1) % Doppler
% figure(1)
% plot((0:length(x(1,:))-1)*Ti,x(3,:),'LineWidth',1) % Der Doppler
% legend('phase','freq','der freq')
% figure(4) % Potencia en Correladores
% hold on
% plot((0:length(Pp)-1)*Ti,abs(Pp),'linewidth',1); grid on; hold on;
% plot((0:length(Ep)-1)*Ti,abs(Ep),'linewidth',1);
% plot((0:length(Lp)-1)*Ti,abs(Lp),'linewidth',1); 
% 
figure(5)
hold on
plot((0:length(Pp)-1)*Ti,real(Pp),'linewidth',1); grid on; 
plot((0:length(Pp)-1)*Ti,imag(Pp),'linewidth',1); 

% figure(6)
% hold on
% plot(CN0_vec(1:end),squeeze(Kk(1,1,1:end-1))') % Ganancia de Kalman para fase
% plot(CN0_vec(1:end),squeeze(Kk(2,1,1:end-1))') % Ganancia de Kalman para frec
% plot(CN0_vec(1:end),squeeze(Kk(3,1,1:end-1))') % Ganancia de Kalman para der frec


% Falta valores para la inicialización
% Falta gráficos y ver como sacamos el exceso de Doppler que es lo
% interesante
% Cambio en el cálculo de la CN0
% Que medidas me llegan y como las puedo poner en las innovations