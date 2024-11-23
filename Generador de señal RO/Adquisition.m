%==========================================================================
%                        Adquisition function GPS
%==========================================================================


function [Y1] = Adquisition(z,Fs,Ti,F_rango,SV,fFI,fig_adq)

% % Parámetros de entrada
% 
% % z señal de entrada
% % Fs frecuencia de muestreo
% % Ti tiempo de adq
% % df paso en frecuencia
% % F_rango frecuencias a recorrer centradas en fFI
% % SV satélite en vista
% % INT_NC cantidad de integraciones no coherentes
% % fig_adq activar figuras 
% 
% % Parámetros de salida 
% 
% % Retardo en chips
% % Frecuencia Doppler
% 
% 
% Ts = 1/Fs; % Tiempo de muestreo
% cx = cacode (SV); % Adquiero un período de 1023 chips para el satélite elegido
% Tchip = 1/(1023e3); % Tiempo de chip nominal 
% C = 3e8; % Velocidad de la luz
% M = floor(Ti/Ts); % Cantidad de muestras a procesar
% K = F_rango/df ; % Cantidad de pasos en frecuencia 
% feini  = fFI-F_rango/2; % Frecuencia de inicio de la busqueda 
% n = 1:M; % Ventana de integración
% c=cx(mod(floor(n*Ts/Tchip),1023)+1); % Replica de código en banda base
% cf=fft(c); %fft de código en banda base
% rzt = 0;
% 
% for i=0:INT_NC-1
%     fe = feini;
%     for k = 1:K
%         zp = z(i*M+1:(i+1)*M).*exp(-1j*(2*pi*fe*n*Ts)); % Quiero correlacionar con una replica
%         % de código local y las señal que llegó
%         zpf = fft(zp);
%         rzs(k,:) = ifft(conj(cf).*zpf); % Función de intercorrelación
%         fe = fe + df;
%     end
%     rzt = rzt + abs(rzs).^2;
% end
% 
% rzs=rzt;
% 
% % Busqueda del máximo 
% [max_value, linear_index] = max(rzs(:));
% [row, col] = ind2sub(size(rzs), linear_index);
% 
% a = (1:M)*Ts/Tchip;
% b= F_rango/2:-df:-F_rango/2;
% 
% 
% 
% retardo_find =a(col+1); % chips de retardo
% frecuencia_find = b(row+1) ; % Frecuencia de pico máximo (no centrada en fFI)
% 
% if fig_adq
% 
%     figure
%     surf((1:M)*Ts/Tchip, -(feini:df:fe-1)+fFI,abs(rzs));
%     shading interp;
%     ylabel('$f-f_0$','Interpreter','latex');
%     xlabel('$\tau [chips]$','Interpreter','latex');
%     zlabel('$|r_{s\widetilde{s}}|$','Interpreter','latex');
%     title('Plano Retardo Doppler','Fontsize',14,'Interpreter','Latex');
% 
% end
% 
% 
% %%

Tchip = 1/(1023e3); % Tiempo de chip nominal 
Ts = 1/Fs;
n = 0:Ti*Fs-1;
paso_Doppler = 0.1*1/Ti;
frecuencias_Doppler = fFI -F_rango/2:paso_Doppler:F_rango/2; % Rango de frecuencias a buscar 
retardos = 0:1/Fs:1e-3-1/Fs; % Un perídodo

N = length(n); % Cantidad de bins en la dimensión del retardo
M = length(frecuencias_Doppler); % Cantidad de bins en la dimensión de frecuencia
K = round(length(z)*Ts/Ti); % Cantidad de integraciones coherentes a realizar

c = cacode(SV);
c = c((mod(floor(n*Ts/Tchip),1023)+1));    
C = conj(fft(c));
Y1=0;
for kk = 1: K

    zq = z(1+kk*N:(kk+1)*N);
    ind_f = 1;

    for ff = frecuencias_Doppler

        y = zq.*exp(-1j*2*pi*ff*n*Ts);
        Y = fft(y);
        r = ifft(Y.*C);
        Y1 = Y1 + abs(r).^2/length(r);
        ind_f = ind_f + 1;
    end
end



if fig_adq
surf(retardos,frecuencias_Doppler,y1); shading interp
end

end