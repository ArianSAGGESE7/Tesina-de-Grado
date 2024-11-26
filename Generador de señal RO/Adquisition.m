%==========================================================================
%                        Adquisition function GPS
%==========================================================================


function [delay,freq,retardos,frecuencias_Doppler,Y1] = Adquisition(z,Fs,Ti,F_rango,SV,fFI,fig_adq,K)



Tchip = 1/(1023e3); % Tiempo de chip nominal 
Ts = 1/Fs;
n = 0:Ti*Fs-1;
paso_Doppler = 0.05*1/Ti;
frecuencias_Doppler = fFI -F_rango/2:paso_Doppler:F_rango/2; % Rango de frecuencias a buscar 
% retardos = 0:1/Fs:Ti-1/Fs; % Un perídodo
retardos = (0:floor(Ti/Ts)-1)*Ts/Tchip;
N = length(n); % Cantidad de bins en la dimensión del retardo
M = length(frecuencias_Doppler); % Cantidad de bins en la dimensión de frecuencia

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
        r(ind_f,:) = ifft(Y.*C);

        ind_f = ind_f + 1;
    end

        Y1 = Y1 + abs(r).^2;
end

% Busqueda del máximo 
[max_value, linear_index] = max(Y1(:));
[row, col] = ind2sub(size(Y1), linear_index);



delay =retardos(col-1); % chips de retardo
freq = -frecuencias_Doppler(row-1) ; % Frecuencia de pico máximo (no centrada en fFI)

if fig_adq
surf(retardos,-frecuencias_Doppler,Y1); shading interp
end

end