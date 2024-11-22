%==========================================================================
%                        Adquisition function GPS
%==========================================================================


function [retardo_find,frecuencia_find] = Adquisition(z,Fs,Ti,df,F_rango,SV,fFI,INT_NC,fig_adq)

% Parámetros de entrada

% z señal de entrada
% Fs frecuencia de muestreo
% Ti tiempo de adq
% df paso en frecuencia
% F_rango frecuencias a recorrer centradas en fFI
% SV satélite en vista
% INT_NC cantidad de integraciones no coherentes
% fig_adq activar figuras 

% Parámetros de salida 

% Retardo en chips
% Frecuencia Doppler


Ts = 1/Fs; % Tiempo de muestreo
cx = cacode (SV); % Adquiero un período de 1023 chips para el satélite elegido
Tchip = 1/(1023e3); % Tiempo de chip nominal 
C = 3e8; % Velocidad de la luz
M = floor(Ti/Ts); % Cantidad de muestras a procesar
K = F_rango/df ; % Cantidad de pasos en frecuencia 
feini  = fFI-F_rango/2; % Frecuencia de inicio de la busqueda 
n = 1:M; % Ventana de integración
c=cx(mod(floor(n*Ts/Tchip),1023)+1); % Replica de código en banda base
cf=fft(c); %fft de código en banda base
rzt = 0;

for i=0:INT_NC-1
    fe = feini;
    for k = 1:K
        zp = 1/sqrt(M)*z(i*M+1:(i+1)*M).*exp(-1j*(2*pi*fe*n*Ts)); % Quiero correlacionar con una replica
        % de código local y las señal que llegó
        zpf = fft(zp);
        rzs(k,:) = ifft(conj(cf).*zpf); % Función de intercorrelación
        fe = fe + df;
    end
    rzt = rzt + 1/INT_NC*abs(rzs).^2;
end

rzs=rzt;

% Busqueda del máximo 
[max_value, linear_index] = max(rzs(:));
[row, col] = ind2sub(size(rzs), linear_index);

a = (1:M)*Ts/Tchip;
b=-(feini:df:fe-df)+fFI;

retardo_find =a(col); % chips de retardo
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