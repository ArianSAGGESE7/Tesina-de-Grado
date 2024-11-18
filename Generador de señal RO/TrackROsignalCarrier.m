%==========================================================================

% Seguimiento de portadora para señal de RO .

%==========================================================================
clc; close all; clear variables

cargardatos = 1;

if cargardatos == 1
    sRO = load("Señal_RO.mat"); % Traemos los datos
    disp(sRO.datosSRO.etiqueta) % Mostramos la metadata correpondiente
end
sROdoppler = sRO.datosSRO.doppler;
sROMuestras = sRO.datosSRO.muestras;
sROamp = sRO.datosSRO.amplitud;
sROfase = sRO.datosSRO.fase;



if sRO.datosSRO.evento == 1
        
    fs = 505001; % Frecuencia de muestreo
    Ts = 1/fs; % Tiempo de muestreo
    tSIM = (0:length(sROMuestras))*Ts; % Tiempo de simulación
    sROdoppler = flip_ro(sRO.datosSRO.doppler);
    sROamp = flip_ro(sRO.datosSRO.amplitud);
    % sROfase = flip_ro(sRO.datosSRO.fase);
    sROMuestras =  flip_ro(sROMuestras);

    % sROMuestras = sROamp.*exp(1j*(2*pi*cumtrapz(0:Ts:tSIM,(sROdoppler))+ sROfase));
    % CN0_db = 45;
    % CN0 = 10^(0.1*CN0_db);
    % N0 = 0.8^2/2/(CN0);
    % N0_var = fs*N0;
    % 
    % % El ruido se genera a partir de una distribución complex normal
    % wI=randn(1,length(sROMuestras));
    % wQ=randn(1,length(sROMuestras));
    % nI=sqrt(N0_var/2).*wI; %Ruido en fase
    % nQ=sqrt(N0_var/2).*wQ; %Ruido en quadratura
    % ruido=nI+1i*nQ; %Ruido "Recibido"
    % 
    % sROMuestras = sROMuestras+ruido;

end



%%

fs = 505001; % Frecuencia de muestreo
Ts = 1/fs; % Tiempo de muestreo
CN0 = 10^(.1*45);
N = length(sROMuestras); % Cantidad de muestras (estos 4 se leen del display)
tSIM = (0:N)*Ts; % Tiempo de simulación
Ti = 5e-3; % Tiempo de integración
N0 = 0.8^2/2/CN0; % N0 fijo

% Parámetros temporales
T = Ti;
M = floor(Ti/Ts) ; % Bloques de muestras de integración 
nt = 0:M-1; % Vecto de puntos que indican cuantas muestras son
MS = floor(N*Ts/Ti); % Cantidad de slots de seguimiento

% Matrices de modelo
F = [1 T T^2/2; 0 1 T; 0 0 1]; % Matriz de transición de estados 
H = [1 T/2 T^2/6]; % Matriz de medición 

% Parametros para las matrices de Kalman --------------------------------------------------------
h0 = 1.8e-20;   % SQGR --> h_0 = 2*10(segs)*(3e-11)^2 con el dato de 3e-11 ADEV @10segs averaging time.      
h_1 = 0;
h2 = 1.24e-21;
q0 = h0/2;
qw = 2*pi^2*h2;
qa = 3e4; % Se obtiene de una tabla que se muestra en libro de Montenbruk;
Q = qa*[T^5/20 T^4/8 T^3/6; T^4/8 T^3/3 T^2/2; T^3/6 T^2/2 T] + qw*[T^3/3 T^2/2 0; T^2/2 T 0; 0 0 0] + q0*[T 0 0; 0 0 0; 0 0 0];
% Q = 0.001*Q;

Q = 0.5*Q;
% % Parametros para las matrices de Kalman --------------------------------------------------------
% h0 = 2e-23;   % SQGR --> h_0 = 2*10(segs)*(3e-11)^2 con el dato de 3e-11 ADEV @10segs averaging time.      
% h_1 = 0;
% h2 = 2.13e-24;
% q0 = h0/2;
% qw = 2*pi^2*h2;
% qa = 3e4; % Se obtiene de una tabla que se muestra en libro de Montenbruk;
% 
% Q = qa*[T^5/20 T^4/8 T^3/6; T^4/8 T^3/3 T^2/2; T^3/6 T^2/2 T]/(3e8)^2 + qw*[T^3/3 T^2/2 0; T^2/2 T 0; 0 0 0]+ q0*[T 0 0; 0 0 0; 0 0 0];
% Q= (2*pi*1545.75e6)^2*Q;
% Q= 0.00001*Q;

% Matriz de Fermin --------------------------------------------------------
% fL1 = 1545.75e6;
% h_1=0;
% h_2=2.13e-24;
% Sj = 90e3;   % Intensidad del jerk. "Smoother-based GPS Signal Tracking in a SDR" pp. 22 da la refe de Sj=1300 rad^2/s^5
% Q = 5e4*[pi^2*T^5/5 pi*T^4/4 pi*T^3/3;pi*T^4/4 T^3/3 T^2/2;pi*T^3/3 T^2/2 T]*Sj;
% % Allan variance covariance model (2-state)
% covxy = (2*pi*fL1)^2*[(h0*T + 2*h_1*T^2 + 2/3*pi^2*h_2*T^3) , (2*h_1*T + pi^2*h_2*T^2) ; ...
%                               (2*h_1*T + pi^2*h_2*T^2)              , (2*h_1 + 2*pi^2*h_2*T)];
% Q = Q + [covxy, [0 0].';[0 0 0]];
% 
% Q = 1*Q;



% Valores que deberian salir de una posterior etapa de tracking
df = 100;
f0 = sROdoppler(1)*(1-0.000009);
% Inicialización de matrices
P0 = [(pi)^2/3 0 0; 0 (2*pi*df)^2/12 0 ; 0 0 100] ; % Matriz de covarianza del error (Hay que incializarla)
x0 = [0;2*pi*(f0); 0]; % Inicialización de los estados con lo que obtuvimos de la etapa de ADQ
x = zeros(3,MS);
Kk = zeros(3,1,MS);    
Px_post = zeros(3,3,MS);
Px_prior = zeros(3,3,MS);
Pp = zeros(1,MS);
Px_post(:,:,1) = P0;
x(:,1)= x0;

% Si bien hay algoritmos que estiman la CN0 cada un determinado tiempo
% nosotros definimos un vector en función de lo conocido

CN0 = sROamp(1:M:length(sROamp)).^2/2/N0;

% Réplica
sRo_replica = exp(1j*(2*pi*f0*nt*Ts)); % Réplica de portadora inicial

for k=0:(MS-2)

    Ssegmento = sROMuestras(M*k +1:M*(k+1)); % Intervalo de muestras 
    
    % Correlador Prompt
    P = sum(conj(sRo_replica).*Ssegmento); %Correlación prompt

    % Discriminador
    Ip = real(P);
    Qp = imag(P);
    D_fase = atan2(Qp,Ip);

    % Estimación por Kalman
    est_x_prior = F*x(:,k+1); % Propagación
    Px_prior(:,:,k+1) = F*Px_post(:,:,k+1)*F.' + Q; % Estimación

    R = 1/(2*Ti*CN0(k+1));
    % R = 1/(2*Ti*CN0);
    K1= H*Px_prior(:,:,k+1)*H.' + R;
    K = Px_prior(:,:,k+1)*H.'/(K1); % Ganancia de Kalman

    est_x_posterior = est_x_prior + K*(D_fase);

    Px_post(:,:,k+2) = Px_prior(:,:,k+1) - K*H*Px_prior(:,:,k+1); % Cálculo de innovaciones


    % Vectores
    x(:,k+2) = est_x_posterior;
    Kk(:,:,k+2) = K; 
    Pp(:,k+1) = P;

    % Generamos réplica local para volver a entrar a los discriminadores 
    xout = F*est_x_posterior;
    sRo_replica = exp(1j*(2*pi*(xout(2)/2/pi)*Ts*nt+xout(1))); % Ajuste de réplica de portadora


end


%% Gráficos
hold on
close all
figure(1)
title('Estados a la salida del filtro','Interpreter','latex')

subplot(3,1,1)
hold on
plot((0:length(Pp)-1)*Ti,atan2(imag(Pp),real(Pp)),'LineWidth',1) % Error de fase
% xlim([0 1])
legend('Salida del discriminador de fase')
subplot(3,1,2)
hold on
plot((0:length(x(1,:))-1)*Ti,x(2,:)/2/pi,'LineWidth',1.5) % Doppler
plot((0:length(sROdoppler)-1)*Ts,sROdoppler,'LineWidth',1) 
% xlim([0 1])
legend('Estimación de Doppler','doppler real')
subplot(3,1,3)
plot((0:length(x(1,:))-1)*Ti,x(3,:)/2/pi,'LineWidth',1) % Doppler-rate
% xlim([0 1])
legend('Doppler rate')


%% Graficos de amplitud y fase del evento 
close all
figure;

subplot(3,1,1)

plot(tSIM(1:end-1),sROamp,'LineWidth',1.5)

subplot(3,1,2)

plot(tSIM(1:end-1),sROfase,'LineWidth',1.5)

subplot(3,1,3)

plot(tSIM(1:end-3),diff(diff(sROfase)./Ts)/Ts,'LineWidth',1.5)

%% Excess Doppler

% El Doppler Geométrico lo interpolamos asi son del mismo largo 

Doppler_Geo_int = interp1((0:length(sROdoppler)-1)*Ts,sROdoppler,(0:length(x(1,:))-1)*Ti,'spline');
figure(2)
plot((0:length(x(1,:))-1)*Ti,Doppler_Geo_int-x(2,:)/2/pi,'LineWidth',1)
title('Excess Doppler')

    % Continuamos en el lazo

    % if k > 8000
    %     c=3e8;
    %     h0 = 2e-23;
    %     h_2=2.13e-24;
    %     qw = 2*pi^2*h_2;
    %     q0 = h0/2;
    % Q = [(T*q0+T^3/3*qw + T^5/20*qa/c^2)  (T^2/2*qw + T^4/8*qa/c^2)  (T^3/6*qa/c^2);
    %  (T^2/2*qw + T^4/8*qa/c^2)  (T*qw + T^3/3*qa/c^2) (T^2/2*qa/c^2)
    %  (T^3/6*qa/c^2)  (T^2/2*qa/c^2) (T*qa/c^2)]*(2*pi*1575.45e6)^2; % Covarianza del proceso de ruido 
    % end