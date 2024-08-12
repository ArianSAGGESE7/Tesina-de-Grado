

function []=Densidad_espectral(x,fs)

    window = hamming(1024*4); % Ventana Hamming de longitud 1024
    noverlap = 512*2; % Mitad de solapamiento
    nfft = 2048*15; % Número de puntos en la FFT
    
    % Aplicar pwelch a la señal s2
    [pxx, f] = pwelch(x, window, noverlap, nfft, fs);
    
    % Ajustar el eje de frecuencias para centrar en fFI
    f_adjusted = f ;%- fs/2;
    
    % Graficar el espectro centrado
    figure;
    plot(f_adjusted, 10*log10(pxx));
    xlabel('Frecuencia (Hz)');
    ylabel('Densidad espectral de potencia (dB/Hz)');
    title('Densidad espectral de potencia de la señal s2 usando pwelch');
    grid on;
    
    % Ajustar los límites del eje x para centrar en fFI
    xlim([-fs/2 fs/2]);
    
end
