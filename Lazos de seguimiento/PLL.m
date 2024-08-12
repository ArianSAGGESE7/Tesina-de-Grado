% Función que realiza las tareas de seguimiento de fase de acuerdo al
% modelo NO LINEAL del PLL.
%
%         [theta] = pll(dt, phi, Kd, Kv, kp, ki)
% dt:     resolución temporal [s]
% phi:    fase de entrada [Hz]
% Kd:     constante del detector de fase [V/rad]
% Kv:     constante del VCO [Hz/V]
% kp, ki: constantes del filtro de lazo F(s) = kp + ki/s;



clear all; clc ; close all
dt = 1e-6;
delta_w=10; % Estcalon de frecuencia
phi =(1:0.01:10)*delta_w; % Rampa de fase de entrada

Kd=4;
Kv=4;
kp=0.19;
ki=.01;

% Valores iniciales de las variables del lazo.
theta = 0;
psi = phi(1) - theta;
ed = Kd*sin(psi);
t = (0:length(phi)-1)*dt; % eje temporal


% Lazo de seguimiento de fase
for tt = 2:length(t)

    psi(tt) = phi(tt-1) - theta(tt-1); % Error de fase

    ed(tt) = Kd*sin(psi(tt)); % Discriminador del error

    e0(tt) = kp*ed(tt) + ki*sum(ed)*dt; % Salida del filtro de lazo

    theta(tt) = Kv*sum(e0)*dt; % La fase es la integral de la salida del filtro

end

figure 
hold on 
grid on
plot(t,ed,'LineWidth',2);
title('salida del discriminador')

figure 
hold on 
grid on
plot(t,e0,'LineWidth',2);
plot(t,delta_w*(ones(1,length(t))),'LineWidth',2);
title('Limite de enganche y salida del filtro de lazo')


figure 
hold on 
grid on
plot(t,theta,'LineWidth',2);
title('fase de salida')


figure 
hold on 
grid on
plot(t,phi,'LineWidth',2);
%end