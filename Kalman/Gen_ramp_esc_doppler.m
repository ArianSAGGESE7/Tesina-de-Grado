
% ========================================================
% Generación de señal con escalones y rampas de frecuencia 
% ========================================================

% Doppler sin aceleración param =0
% Si se quiere escalon de freq param =1
% Si se quiere rampa de freq param =2
% t1 comienzo del evento
% t2 fin del evento 
% Amp amplitud del escalon o amplitud final de la rampa (se mide en aceleraciones g)
function [doppler] = Gen_ramp_esc_doppler(doppler,TD,Ts,param,t1,t2,Amp,Amp_ini)

% doppler value
% TD simulation time
% Ts sample time

n = 0:TD/Ts-1; % Indice de largo simulación
doppler = doppler*ones(1,length(n));

Doppler_end = Amp; % Valor final del Doppler relativo al Doppler definido

if param ==1
    doppler(t1/Ts:t2/Ts) = doppler(1)*Doppler_end; % escalón

elseif param == 2
    doppler(t1/Ts:t2/Ts)= doppler(1)*Amp_ini + doppler(1)*(Doppler_end-1).*(0:length(doppler(t1/Ts:t2/Ts))-1)/(length(doppler(t1/Ts:t2/Ts))-1); % Rampa

  
end