%Lazo Wide Operativo 10s
clear all; close all; clc;
%%
%PANEL DE CONTROL --------------------------------------------------------%
CARGAR_SENAL_SINTETICA = 1; %Habilitar creación y carga de señ. sintética % 
CARGAR_SENAL_REAL      = 0; %Habilitar carga de señal real ---------------%
ADQUIRIR               = 1; %Habilitar adquisición -----------------------%       
SEGUIMIENTO            = 1; %Habilitar seguimiento -----------------------%
    FC_ASIST           = 1; %Habilita asistencia e/ lazos ----------------%              
    CODE_LOOP          = 1; %Habilita lazo de código ---------------------%             
    FREQ_LOOP          = 1; %Habilita lazo de portadora ------------------%
TCH_ADULT              = 1; %Habilita adulteración de tiempo de chip -----%
MC                     = 0; %Habilita rayo de multicamino ----------------%
    taum   =  0.2; %Retardo relativo del rayo de multicamino [chips]------%
    am_1   =  0.5; %Amplitud relativa del rayo de multicamino ------------%
   
%Sobre la señal sintética ------------------------------------------------%
doppler  = 100; %Doppler deseado (en Hz) ---------------------------------%                           
CN0dB    = 45;  %Seleccione relación señal a ruido en dB -----------------%                      
                      
%Código común    
NUMERO_DE_SATELITE = 5;  %Seleccione sat. SV: 16(m),22(m),03(m),23,07,09,06           
cx=cacode(NUMERO_DE_SATELITE); %------------------------------------------%
fL1 = 1575.42e6; %Frecuencia nominal -------------------------------------%
fOL = 1590e6; %Frecuencia de oscilador local -----------------------------%             
fFI = fOL-fL1; %Frecuencia intermedia ------------------------------------%       
Tch = 1e-3/1023; %Tiempo de chip nominal ---------------------------------%
C   = 3e8; %Velocidad de la luz ------------------------------------------%
lambda = C/fL1; %Longitud de onda nominal --------------------------------%
fs = 40e6;  %Frecuencia de muestreo del proyecto --------------------------%
Ts = 1/fs; %Tiempo de muestreo del proyecto ------------------------------%
%%
if CARGAR_SENAL_REAL==1
fileID = fopen('GPSL1-6MHz-14-55.txt');
data   = textscan(fileID,'%f,%f','HeaderLines',1);
fclose(fileID);
%time = cell2mat(data(1));
x    = cell2mat(data(2));
% D0: I0 es cable marron   (SIGNO) -> ch0 MSB. /// SIGNO = 0 (positivo)
% D1: I1 es cable rojo     (MAGN)  -> ch1 LSB. /// MAGN  = 0 (es 1)
% D2: Q0 es cable amarillo (SIGNO) -> ch0 MSB. /// SIGNO = 0 (positivo)
% D3: Q1 es cable naranja  (MAGN)  -> ch1 LSB. /// MAGN  = 0 (es 1)
valores   = [1 3 -1 -3];
Iaux = mod(x,4) + 1;
I = valores (Iaux);
Qaux = bitshift(mod(x,16),-2) + 1;
Q = valores (Qaux);
z=I+1j*Q;
TD=length(z)*Ts; %Duración de datos [seg]
N=sum(abs(z).^2)/length(z); %Estimación de potencia de ruido
end

if CARGAR_SENAL_SINTETICA==1
fdata=50; %Tasa de datos
Tdata=1/fdata; %Tiempo de bit de datos
TD=2; %Seleccione duración de datos [seg]                             
n=0:TD/Ts-1; %Indice largo
PEND=-doppler*lambda*Ts; %Cambio muestra a muestra                                                               
x=500+[1:length(n)]*PEND; %Pseudorango (en metros)
taut=x/C; %Tiempo asociado a pseudorango
cs=cx(mod(floor((n*Ts-taut)/Tch),length(cx))+1); 
ndata=0:TD/Tdata-1; %Indice para datos
data=sign(rand(1,length(ndata))-.5);
cdata=data(mod(floor((n*Ts-taut)/Tdata),length(data))+1);
s1=cdata.*cs.*exp(1j*(2*pi*(-fL1)*(Ts*n-taut)));
s2=s1.*exp(-1j*(2*pi*fOL*Ts*n)); %Demodulación a frecuencia intermedia

%Rayo de multicamino
cdatam=data(mod(floor((n*Ts-taut-taum*Tch)/Tdata),length(data))+1);
csm=cx(mod(floor((n*Ts-taut-taum*Tch)/Tch),length(cx))+1); 
s1m=cdatam.*csm.*exp(1j*(2*pi*(-fL1)*(Ts*n-taut-taum*Tch)));
s2m=am_1*s1m.*exp(-1j*(2*pi*fOL*Ts*n)); %Demodulación a FI

%Sobre el ruido
CN0=10^(.1*CN0dB);                                
N=1/(Ts*CN0);
wI=randn(1,length(n));
wQ=randn(1,length(n));
nI=sqrt(N/2).*wI;   %Ruido en fase
nQ=sqrt(N/2).*wQ;   %Ruido en quadratura
ruido=nI+1i*nQ;  %Ruido "Recibido"

z=s2+MC*s2m+ruido; %Señal total "recibida"
end
%FIN CARGA DATOS-----------------------------------------------------------

%% ADQUISICIÓN---------------------------------------------------------------

TD=length(z)*Ts; %Duración de datos [seg]
N=sum(abs(z).^2)/length(z); %Estimación de potencia de ruido

if ADQUIRIR==1
M=fs/500; %Cantidad de muestras a procesar
T=M*Ts; %Tiempo de integración
F=10e3; %Rango de frecuencia a recorrer
df=0.1*(1/T); %Paso en frecuencia
f_central=fFI; %frecuencia central de búsqueda
K=F/df; %Pasos en frecuencia
feini=f_central-F/2; %Frecuencia inicial

Tch_=Tch;
if TCH_ADULT==1
Tch_=1/(1/Tch-1/(2*T)); %Tiempo de chip adulterado
end

n=0:M-1; %Indice para ventana de integración
c=cx(mod(floor(n*Ts/Tch_),1023)+1); %Código en banda base
cf=fft(c); %fft de código en banda base

zp1=z(1:M);  %Primer intervalo de datos
zp2=z(M+1:2*M); %Segundo intervalo de datos
zp3=z(2*M+1:3*M); %Tercer intervalo de datos

UMB_co=5*sqrt(M*N); %Umbral de correlación 

%Inicio Plano 1
fe=feini;
for k=1:K  %Barrido en frecuencia
zpb=zp1.*exp(-1j*(2*pi*fe*n*Ts)); %Datos demodulados
zpbf=fft(zpb); %fft de los datos en banda base      
rzs(k,:)=ifft(conj(cf).*zpbf); %Función de intercorrelación entre la parte de señal recibida y la réplica propuesta
fe=fe+df;
end
rzs1=abs(rzs).^2;
%Fin Plano 1

%Inicio Plano 2
fe=feini;
for k=1:K  %Barrido en frecuencia
zpb=zp2.*exp(-1j*(2*pi*fe*Ts*n)); %Datos demodulados
zpbf=fft(zpb); %fft de los datos en banda base      
rzs(k,:)=ifft(conj(cf).*zpbf); %Función de intercorrelación entre la parte de señal recibida y la réplica propuesta
fe=fe+df;
end
rzs2=abs(rzs).^2;
%Fin Plano 2

%Inicio Plano 3
fe=feini;
for k=1:K  %Barrido en frecuencia
zpb=zp3.*exp(-1j*(2*pi*fe*Ts*n)); %Datos demodulados
zpbf=fft(zpb); %fft de los datos en banda base      
rzs(k,:)=ifft(conj(cf).*zpbf); %Función de intercorrelación entre la parte de señal recibida y la réplica propuesta
fe=fe+df;
end
rzs3=abs(rzs).^2;
%Fin Plano 3

%Inicio Plano Total
rzs_total=rzs1+rzs2+rzs3; %Plano total
fe=feini;
PC=0; %Pico de correlación inicializado
for k=1:K  %Barrido en frecuencia
if(max(abs(rzs_total(k,:)))>PC)  %Si el máximo de correlación actual supera al candidato anterior
    PC=max(abs(rzs_total(k,:)));  %Actualizo con nuevo candidato
    fxe=fe;  %Actualizo frecuencia candidata
    tau=find(abs(rzs_total(k,:))==PC)-1;  %Candidatos para retardo de código, en chips
    tau_=tau(1)*(tau(1)<6000)+(tau(1)-6000)*(tau(1)>=6000);
    tauxe=(tau_-1)*Ts/Tch;  %Me quedo con el primero
end
fe=fe+df;
end

if(PC>UMB_co) %Si el candidato encontrado supera el umbral
    fx_prev=fxe; %Declaro frecuencia
    taux_prev=tauxe; %Declaro retardo (en chips)     
end
%Fin primera adquisición. En fx_prev, taux_prev está lo adquirido

%Inicia refinamiento de frecuencia
M=fs/100; %Cantidad de datos a procesar
T=M*Ts; %Tiempo de integración
nt=0:M-1;
F=200; %Rango de frecuencia a recorrer
df_fino=0.01*(1/T); %Paso en frecuencia
K=F/df_fino; %Cantidad de pasos en frecuencia
feini=fx_prev-F/2; %Frecuencia de escaneo inicial
dt_fino=.5; %Paso en retardo (en chips)
nc=11; %Número de chips a recorrer (debe ser impar)
tauxini=taux_prev-(nc-1)/2;
zp=z(1:M);  %Selecciono intervalos de datos

Tch_=Tch;
if TCH_ADULT==1
Tch_=1/(1/Tch-1/(2*T)); %Tiempo de chip adulterado
end

%Creación de plano ultra fino
fe=feini;
PC=0;
for k=1:K
    s=exp(1j*(2*pi*fe*nt*Ts)); %Réplica de portadora 
    taux=tauxini;
    for i=1:nc/dt_fino   
    c=cx(mod(floor((nt*Ts-taux*Tch_)/Tch_),1023)+1); %Réplica de código
    rzs_fino(k,i)=sum(conj(c.*s).*zp);
    if abs(rzs_fino(k,i))>PC
        PC=abs(rzs_fino(k,i));
        fxe=fe;
        tauxe=taux;
    end
    taux=taux+dt_fino;
    end
    fe=fe+df_fino;
end

UMB_co=2*sqrt(M*N); %Umbral de correlación

if(PC>UMB_co) %Si el candidato encontrado supera el umbral
    ADQ_OK=1;
    fx_fino=fxe; %Declaro frecuencia
    taux_fino=tauxe; %Declaro retardo (en chips)
    x_fino=taux_fino*1540*C/fL1;  %Pseudorango estimado por la adquisición (opcional)
end
%Fin refinamiento de frecuencia
end
%FIN ADQUISICIÓN-----------------------------------------------------------

%% BUSQUEDA DE BIT DE DATOS--------------------------------------------------
M=fs/100; %Cantidad de datos a procesar
T=M*Ts; %Tiempo de integración
ntt=0:M-1;
c=cx(mod(floor((ntt*Ts-taux_fino*Tch_)/Tch_),length(cx))+1); %Réplica de código
s=exp(1j*(2*pi*fx_fino*ntt*Ts)); %Réplica de portadora 
MSS=(TD/1e-3)-10;  %Cantidad de correlaciones a hacer 
CDM=1e-3/Ts;  %Cantidad de datos por milisegundo (creyéndole al reloj)
fdoppler=fFI-fx_fino;  
delta_T_frec=1/(fdoppler+fL1)-(1/fL1); %Delta de cambio del periodo de portadora
delta_T_code=1023*1540*delta_T_frec;  %Retardo o deriva de codigo debido al doppler 

for k=0:MSS-1
zvd=z(CDM*k+1:CDM*(k+10)); %Selección de slot        
P_VD(k+1)=abs(sum(conj(c.*s).*zvd)); %Correlación    
c=cx(mod(floor((ntt*Ts-taux_fino*Tch_-(k+1)*delta_T_code)/Tch_),length(cx))+1); %Réplica de código
end

PROMPE=0.5*sum(P_VD)/length(P_VD);
plot(0:MSS-1,P_VD)
for k=1:length(P_VD)
ANTERIOR=P_VD(k);
CENTRAL=P_VD(k+1);
POSTERIOR=P_VD(k+2);
if ((CENTRAL<POSTERIOR)&&(CENTRAL>ANTERIOR)&&(CENTRAL<PROMPE))
MOM_TRAN=(k+4)*1e-3; %Instante en que ocurre la primera transición de bit de datos, tomando como referencia el origen de la ristra generada
break
end
end
%FIN BUSQUEDA DE BIT DE DATOS----------------------------------------------

%% SEGUIMIENTO---------------------------------------------------------------
fx_fino=fFI+100;
if SEGUIMIENTO==1

    fx=fx_fino; %Inicializo con frecuencia adquirida [Hz]
    taux=taux_fino+(fx-fFI)/1540*MOM_TRAN; %Inicializo con retardo de código adquirido [chips]

    %Para el trackeo:
    M=fs/100; %Cantidad de datos a procesar por cada ciclo de seguimiento
    T=M*Ts; %Tiempo de integración
    nt=0:M-1;
    MS=floor((TD-MOM_TRAN)/T);  %Cantidad de slots de seguimiento

    Tch_=Tch;
    if TCH_ADULT==1
        Tch_=1/(1/Tch-1/(2*T)); %Tiempo de chip adulterado
    end

    %Dinámica de lazos
    BW_code=1; %BW de filtro de lazo de código                              <<
    K0=4*BW_code*T; %Cte de lazo de código

    BW_freq=15;  %BW de filtro de lazo de portadora                          <<
    a2=2.4*BW_freq/0.7845*T; %Ctes de lazo de portadora
    a1=1.1*(BW_freq/0.7845*T)^2;
    a0=(BW_freq/0.7845*T)^3;

    out_p1=0;
    out_p0=(fx-fFI)*2*pi*T;

    c=cx(mod(floor((nt*Ts-taux*Tch_)/Tch_),length(cx))+1); %Réplica de código inicial
    s=exp(1j*(2*pi*fx*nt*Ts)); %Réplica de portadora inicial

    %Inicio de trackeo
    for k=0:MS-1
        zseg=z(M*k+MOM_TRAN/Ts:M*(k+1)+MOM_TRAN/Ts-1); %Selección de slot

        %Lazo de portadora:
        P=sum(conj(c.*s).*zseg); %Correlación prompt
        IPS=real(P);
        QPS=imag(P);
        D=atan(QPS/IPS); %Discriminador
        out_p0=out_p0+a2*D; %variable de estado
        out_p1=out_p1+a1*D+out_p0; %variable de estado
        fi_e=a0*D+out_p1; %Salida del lazo: fase estimada
        fdop=out_p0/(2*pi*T); %Frecuencia estimada

        %Lazo de código:
        sE=s.*cx(mod(floor((nt*Ts-CODE_LOOP*taux*Tch_+.5*Tch_)/Tch_),1023)+1); %Réplica early
        E=sum(conj(sE).*zseg); %Correlación early
        sL=s.*cx(mod(floor((nt*Ts-CODE_LOOP*taux*Tch_-.5*Tch_)/Tch_),1023)+1); %Réplica late
        L=sum(conj(sL).*zseg); %Correlación late
        delta_tau=.5*(abs(E)-abs(L))/(abs(E)+abs(L)); %Discriminador
        f_tau=(fdop/1540); %Frecuencia del código alterada por doppler [chips/seg]

        taux=taux-K0*delta_tau+FC_ASIST*f_tau*T;  %Variable de estado del lazo de código y asistencia (realimentación)
        %del lazo de portadora
        s=exp(1j*(2*pi*(fFI+FREQ_LOOP*fdop)*Ts*nt+FREQ_LOOP*fi_e)); %Ajuste de réplica de portadora
        c=cx(mod(floor((nt*Ts-CODE_LOOP*taux*Tch_)/Tch_),1023)+1); %Ajuste de réplica de código

        x_seg(k+1)=taux*Tch*C; %Estimación de pseudorango
        fDD(k+1)=-fdop; %Estimación de doppler
        pE(k+1)=E; %Potencia de correlación early
        pL(k+1)=L; %Potencia de correlación late
        pP(k+1)=P; %Potencia de correlación prompt
    end
end
%%
close all
% figure
% %pseudorango
% subplot(2,2,1); plot(1:10:10*length(x_seg),x_seg,'linewidth',1); grid on; 
% title(['Pseudorango estimado - SV',num2str(NUMERO_DE_SATELITE)],'Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
% ylabel('Pseudorango [m]','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
% xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% 
% %frecuencia
% subplot(2,2,2); plot(1:10:10*length(fDD),fDD,'linewidth',1); grid on;
% title(['Doppler estimado - SV',num2str(NUMERO_DE_SATELITE)],'Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% ylabel('Frecuencia [Hz]','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
% xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
figure(1)
%Potencias de correlación
plot(1:10:10*length(pP),abs(pP),'linewidth',1); grid on; hold on;
plot(1:10:10*length(pE),abs(pE),'linewidth',1); hold on
plot(1:10:10*length(pL),abs(pL),'linewidth',1); 
title(['Potencias P-E-L - SV',num2str(NUMERO_DE_SATELITE)],'Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
ylabel('Potencia [-]','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
figure(2)
%I y Q de potencia prompt
plot(1:10:10*length(pP),abs(real(pP)),'linewidth',1); grid on; hold on
plot(1:10:10*length(pP),imag(pP),'linewidth',1);
title(['Potencias I-Q Prompt - SV',num2str(NUMERO_DE_SATELITE)],'Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
ylabel('Potencia [-]','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');

% figure
% plot(1:length(P_VD),P_VD); grid on
% title('Correlaciones por ventana deslizante','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% ylabel('|r(s,z)|','Fontsize',14,'FontAngle','italic','Interpreter','Latex'); 
% xlabel('Tiempo [ms]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% 
% figure
% surf(abs(rzs_total));
% shading interp; 
% ylabel('Proporcional al Doppler','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% xlabel('Proporcional al retardo de codigo [chips]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% zlabel('|Rzs|','Fontsize',14,'FontAngle','italic','Interpreter','Latex');   
% title(['Plano de primera Adquisicion - SV',num2str(NUMERO_DE_SATELITE)],'Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% 
% figure
% surf(abs(rzs_fino));
% shading interp; 
% ylabel('Proporcional al Doppler','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% xlabel('Proporcional al retardo de codigo [chips]','Fontsize',14,'FontAngle','italic','Interpreter','Latex');
% zlabel('|Rzs|','Fontsize',14,'FontAngle','italic','Interpreter','Latex');   
% title(['Plano del Refinamiento - SV',num2str(NUMERO_DE_SATELITE)],'Fontsize',14,'FontAngle','italic','Interpreter','Latex');