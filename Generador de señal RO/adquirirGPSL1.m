function [resultados_L1] = adquirirGPSL1(datos,fs,fol1,dur,graficar)



fc = 1.023e6; % Frecuencia de chip
dwell = 1e-3; % Tiempo de integración coherente para L1
Pfa_target_L1 = 1e-3; % Probabilidad de falsa alarma objetivo para L1

satL1_buscar = 3; % Vector que contiene la lista de satélites a buscar en L1


z1 = datos(:)';

FREQ_GPS_L1 = 1.57542e9;

f_FI = FREQ_GPS_L1 - fol1;

n = 0:dwell*fs-1;

paso_dopp = 0.1*1/dwell; % Paso de búsqueda en la dimensión Doppler
dopp_max = 40e3; % Máximo Doppler a buscar
dopp_central = 0; % Frecuencia central de desviación Doppler
frecs_dopp = dopp_central-dopp_max:paso_dopp:dopp_central+dopp_max; % Frecuencias Doppler a buscar
rets = 0:1/fs:1e-3-1/fs; % Retardos a buscar (para tener el dato en seg.)

N = length(n); % Cantidad de bins en la dimensi�n retardo
M = length(frecs_dopp); % Cantidad de bins en la dimensión frecuencia con rango -10kHz-10kHz
K = round(dur/dwell); % Cantidad de integraciones coherentes
P = N*M*K-K;

x = (1:10000)*.005*(K/P)*(N*M);

pfa = 1-(fcdf(x,2*K,2*P)).^(N*M); % Pfa para detector cuadr�tico full parallel CFAR

thresh = min(x(pfa<Pfa_target_L1)); % Umbral de detección para la Pfa_target

satsL1_vista = [];
doppler_L1 = [];
retardo_L1 = [];
CN0_L1 = [];

amp_OL = sqrt(2);
y1=0;
for svnum = satL1_buscar
        
    c = cacode(svnum);
    c = c((mod(floor(n*fc/fs),1023)+1));    
    C = conj(fft(c));
    
    
    for kk = 0:3
        zq = z1(1+kk*N:(kk+1)*N);
%         z_dem = zq.*exp(1i*2*pi*(f_FI)*n/fs);
        
        I1 = zq.*amp_OL.*cos(2*pi*(f_FI)*n/fs);
        Q1 = zq.*amp_OL.*sin(2*pi*(f_FI)*n/fs);

        ind_f = 1;
        for ff = frecs_dopp
            y = (I1+1i*Q1).*exp(-1i*2*pi*ff*n/fs);
            Y = fft(y);
            r(ind_f,:) = ifft(Y.*C);
            ind_f = ind_f+1;
        end

        y1 = y1 + abs(r/256).^2/length(r); 
    end
    

    y1 = y1(:,1:fs*1e-3);
    [m1, i1] = max(y1(:));
    media1 = mean(y1(:));
    [frec1, ret1] = ind2sub(size(y1),i1);
    CN0_L1_est = 10*log10(1/(2*dwell)*m1/mean(y1(:)));
    
    if (m1/media1)>thresh
        satsL1_vista = [satsL1_vista svnum];
        doppler_L1 = [doppler_L1 frecs_dopp(frec1)];
        retardo_L1 = [retardo_L1 rets(ret1)];
        CN0_L1 = [CN0_L1 CN0_L1_est];
        % Para graficar solamente los visibles
        if graficar
            figure('units','normalized','outerposition',[0 0 1 1]) % Creo y maximizo figura.
            surf(rets,frecs_dopp,y1);title("Satélite GPS Nº" + svnum + " con Integracion de " + dur*1000 + " ms"); shading interp
            xlabel('Retardo [s]'); ylabel('Desviaci�n Doppler [Hz]');            
%             figure; plot(frecs_dopp,y1(:,ret1))
            pause
%             close;
        end
    end
    
end


resultados_L1 = [satsL1_vista' doppler_L1' retardo_L1' CN0_L1'];

satsL1_vista
CN0_L1



% figure; plot(frecs_dopp,y1(:,ret1))
end