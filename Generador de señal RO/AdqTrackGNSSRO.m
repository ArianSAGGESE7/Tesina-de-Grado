%--------------------------------------------------------------------------
%        RADIO OCCULTATION GNSS - TRACKING METHOD USING KALMAN FILTER    
%--------------------------------------------------------------------------

% Los datos para correr este script se obtienen de GenParamRO
clc; close all; clear variables
cargardatos = 1;

if cargardatos == 1
    sRO = load("Señal_RO.mat"); % Traemos los datos
    disp(sRO.datosSRO.etiqueta) % Mostramos la metadata correpondiente
    % SigRoMuestras = readVector('SigRoGnss');
end
SigRoMuestras = sRO.datosSRO.muestras;
sROdoppler = sRO.datosSRO.doppler;
sROamp = sRO.datosSRO.amplitud;
sROfase = sRO.datosSRO.fase;
fs = sRO.datosSRO.fs;
CN0db = sRO.datosSRO.CN0 ;

Ts = 1/fs;
Tin = 0;
Tsim = length(SigRoMuestras)*Ts;
tInterp= (Tin:Ts:Tsim); % Tiempo a evaluar (con la resolución que se requiera, esto lo fijamos nosotros)
SV = 1; % Adquisición del satélite 
cx = cacode (SV); % Adquiero un período de 1023 chips para el satélite elegido
fL1 = 1575.42e6; % Frecuencia nominal GPS
fOL = 1575.42e6; % Frecuencia de oscilador local (no es necesariamente la misma)
fFI = fL1-fOL; % Frecuencia intermedia
Tchip = 1/(1023e3); % Tiempo de chip nominal 
C = 3e8; % Velocidad de la luz


%% Etapa de Adquisición

Tadq = 1e-3; % Tiempo de integración de adquisición
F_rango = 90e3; % Rango de frecuencias a recorrer
df = 0.05/Tadq; % Este paso en frecuencia tiene que ver con el error máximo que toleramos en el Kalman
INT_NC = 3; % Integraciones no coherentes
fig_adq =1; % Habilitar figuras
K = 1; % Promediaciones

[ret0,frec0] = Adquisition(SigRoMuestras,fs,Tadq,F_rango,SV,fFI,fig_adq,K)

%% Perdida de potencia por no actualizar el retardo

Tinte = 10e-3; % Tiempo de integración
ntt = 0:Tinte/Ts-1; %Indice para ventana de integración de T segundos
MSS = floor(length(SigRoMuestras)*Ts/Tinte); %Cantidad de correlaciones a hacer
CDM = floor(Tinte/Ts); %Cantidad de datos por milisegundo a una tasa de muestreo de fs, entran CMD datos.
c=cx(mod(floor((ntt*Ts-ret0*Tchip)/Tchip),length(cx))+1);% Replicas de portadora y de código
s=exp(-1j*(2*pi*frec0*ntt*Ts)); % no cambia la frecuencia por lo que se actualiza una única vez.
fdoppler=frec0-fFI;
delta_T_frec=1/(fdoppler+fL1)-(1/fL1);% Ayuda del lazo de portadora al de código, si no lo tenemos perfemos% potencia
delta_T_code=1023*1540*delta_T_frec;
 P_VD = zeros(1,MSS);
for k=0:MSS-1
 zvd=SigRoMuestras(CDM*k+1:CDM*(k+1)); %Selección de slot
 P_VD(k+1)=abs(sum(conj(c.*s).*zvd)); %Correlación
 c=cx(mod(floor((ntt*Ts-ret0*Tchip-k*delta_T_code)/Tchip),length(cx))+1); %Réplica de código
end

plot((1:length(P_VD))*Tinte,P_VD);



%% Etapa de Seguimiento 

TranBit = 5E-3; % corrimiento dentro del slot para no tomar un bit de datos a la mitad
Ti = 1e-3; % Tiempo de integración

CN0 = 10^(.1*CN0db);
N = length(SigRoMuestras); % Cantidad de muestras (estos 4 se leen del display)
tSIM = (0:N)*Ts; % Tiempo de simulación
N0 = 1^2/2/CN0; % N0 fijo
M = floor(Ti/Ts) ; % Bloques de muestras de integración 
nt = 0:M-1; % Vecto de puntos que indican cuantas muestras son
MS = floor((N*Ts-TranBit)/Ti); % Cantidad de slots de seguimiento
tKalman = (1:MS)*Ti;

% Matrices de Kalman
F = [1 Ti Ti^2/2; 0 1 Ti; 0 0 1]; % Matriz de transición de estados 
H = [1 Ti/2 Ti^2/6]; % Matriz de medición 

% Parametros para las matrices de Kalman --------------------------------------------------------
h0 = 1.8e-20;   % SQGR --> h_0 = 2*10(segs)*(3e-11)^2 con el dato de 3e-11 ADEV @10segs averaging time.      
h_1 = 0;
h2 = 1.24e-21;
q0 = h0/2;
qw = 2*pi^2*h2;
qa = 3e4; % Se obtiene de una tabla que se muestra en libro de Montenbruk;
Q = qa*[Ti^5/20 Ti^4/8 Ti^3/6; Ti^4/8 Ti^3/3 Ti^2/2; Ti^3/6 Ti^2/2 Ti] + qw*[Ti^3/3 Ti^2/2 0; Ti^2/2 Ti 0; 0 0 0] + q0*[Ti 0 0; 0 0 0; 0 0 0];
% Q = 0.001*Q;

Q = .5*Q;
% Inicializamos Kalman y lazo de código.
f0 = frec0;

taux = ret0 - f0/1540*TranBit; %==

% Parámetros del filtro de lazo
BW_code = 20; % Ancho de banda del filtro en [Hz]
K0 = 4*BW_code*Ti; % Constante de lazo de código

% Inicialización de matrices
P0 = [(pi)^2/3 0 0; 0 (2*pi*df)^2/12 0 ; 0 0 100] ; % Matriz de covarianza del error (Hay que incializarla)
x0 = [0;2*pi*(f0); 0]; % Inicialización de los estados con lo que obtuvimos de la etapa de ADQ
x = zeros(3,MS);
Kk = zeros(3,1,MS);    
Px_post = zeros(3,3,MS);
Px_prior = zeros(3,3,MS);
Pp = zeros(1,MS);
Ep = zeros(1,MS);
Lp = zeros(1,MS);
Px_post(:,:,1) = P0;
x(:,1)= x0;

% CN0 = sROamp(1:M:length(sROamp)).^2/2/N0; % Definida y tomada como conocida
CN0 = ones(1,MS)*10^(0.1*CN0db);

c = cx(mod(floor((nt*Ts-taux*Tchip)/Tchip),length(cx))+1); %Réplica de código inicial
s = exp(-1j*(2*pi*f0*nt*Ts)); %Réplica de portadora inicial


for k=0:(MS-2)

    zseg = SigRoMuestras(M*k + floor(TranBit/Ts) :M*(k+1) + floor(TranBit/Ts)-1); % Intervalo de muestras 
    
    % Correladores Prompt, Early y Late
    P = sum(conj(c.*s).*zseg); %Correlación prompt
    sE = s.*cx(mod(floor((nt*Ts-taux*Tchip+.5*Tchip)/Tchip),1023)+1); %Réplica early
    E = sum(conj(sE).*zseg); %Correlación early
    sL = s.*cx(mod(floor((nt*Ts-taux*Tchip-.5*Tchip)/Tchip),1023)+1); %Réplica late
    L=sum(conj(sL).*zseg); %Correlación late

                                % Discriminadores
    Ip = real(P);
    Qp = imag(P);
    %Discriminador de fase
    D_fase = atan(Qp/Ip);
    %Discriminador de código
    D_tau=.5*(abs(E)-abs(L))/(abs(E)+abs(L));
    
    

    % Estimación por Kalman
    est_x_prior = F*x(:,k+1); % Propagación
    Px_prior(:,:,k+1) = F*Px_post(:,:,k+1)*F.' + Q; % Estimación

    R = 1/(2*Ti*CN0(k+1));
  
    K1= H*Px_prior(:,:,k+1)*H.' + R;
    K = Px_prior(:,:,k+1)*H.'/(K1); % Ganancia de Kalman

    est_x_posterior = est_x_prior + K*(D_fase);

    Px_post(:,:,k+2) = Px_prior(:,:,k+1) - K*H*Px_prior(:,:,k+1); % Cálculo de innovaciones


    % Vectores
    x(:,k+2) = est_x_posterior;
    Kk(:,:,k+2) = K; 
    Pp(:,k+1)=P;
    Ep(:,k+1)=E;
    Lp(:,k+1)=L;
    
    xout = F*est_x_posterior; % Propagamos una vez mas
    %Estimación del retardo
    taux= taux - K0*D_tau + xout(2)*Ti/1540/2/pi + xout(3)*Ti^2/2/1540/2/pi; 
    
    % Generamos réplica local para volver a entrar a los discriminadores 
    s = exp(-1j*(2*pi*(fFI+xout(2)/2/pi)*Ts*nt + xout(1))); % Ajuste de réplica de portadora
    c = cx(mod(floor((nt*Ts-taux*Tchip)/Tchip),1023)+1); % Ajuste de réplica de código

    % Continuamos en el lazo
end


% Gráfico (Estados, dicriminadores y correladores)
% close all
hold on

figure(1)
title('Estados a la salida del filtro','Interpreter','latex')
subplot(3,1,1)
hold on
plot((1:length(Pp))*Ti,atan(imag(Pp)./real(Pp)),'LineWidth',1) % Error de fase
legend('Salida del discriminador de fase')
subplot(3,1,2)
hold on
plot((0:length(x(1,:))-1)*Ti,x(2,:)/2/pi,'LineWidth',1) % Doppler
plot((0:length(sROdoppler)-1)*Ts,sROdoppler)
legend('Estimación de Doppler','Doppler real')
subplot(3,1,3)
plot((0:length(x(1,:))-1)*Ti,x(3,:)/2/pi,'LineWidth',1) % Doppler-rate
legend('Doppler rate')

figure(2) 
hold on
plot((1:length(Pp))*Ti,abs(Pp),'linewidth',1); grid on;
plot((1:length(Pp))*Ti,abs(Ep),'linewidth',1);
plot((1:length(Pp))*Ti,abs(Lp),'linewidth',1); 
legend('Prompt','early','late','Interpreter','latex')

% figure(3)
% Doppler_Geo_int = interp1((0:length(sROdoppler)-1)*Ts,sROdoppler,(0:length(x(1,:))-1)*Ti,'spline');
% plot((0:length(x(1,:))-1)*Ti,medfilt1(Doppler_Geo_int-x(2,:)/2/pi,50),'LineWidth',1)
% title('Excess Doppler')