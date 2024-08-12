close all;
clear all;
% Registros de acumulación para integración
reg1=0;
reg2=0;
reg3=0;
eta=sqrt(2)/2; % Contante de amortiguación PLL 
theta=2*pi*1/100; % Phase
kp=(4*eta*theta)/(1+2*eta*theta+theta^2); % Constante de diseño del filtro de lazo
ki=(4*theta^2)/(1+2*eta*theta+theta^2); % Constante de diseño del filtro de lazo
d_phi_1=1/50; % Paso de fase para la simulación
n_data=1000;
for nn=1:n_data
    phi1=reg1+d_phi_1; % Incrementamos la fase para la simulación
    phi1_reg(nn)=phi1;
    s1=exp(1i*2*pi*reg1); % Señal de entrada
    s2=exp(1i*2*pi*reg2); % Réplica local
    s1_reg(nn)=s1;
    s2_reg(nn)=s2;
    t=s1*conj(s2); % Producto que entra al discriminador de fase
    phi_error=atan(imag(t)/real(t))/(2*pi); % Discriminador
    phi_error_reg(nn)=phi_error;
    sum1=kp*phi_error+phi_error*ki+reg3; % Integral de los filtros
    reg1_reg(nn)=reg1;
    reg2_reg(nn)=reg2;
    reg1=phi1;
    reg2=reg2+sum1;
    reg3=reg3+phi_error*ki;
    phi2_reg(nn)=reg2;
end
figure(1)
plot(phi1_reg);
hold on
plot(phi2_reg,'r')
hold off;
grid on;
title('Fase');
xlabel('Muestras');
ylabel('Fase');
figure(2)
plot(phi_error_reg);
title('Erorr del detector de fase');
grid on;
xlabel('muestras(n)');
ylabel('Error de fase(Grados)');
figure(3)
plot(real(s1_reg));
hold on;
plot(real(s2_reg),'r');
hold off;
grid on;
title('Señal de entrada y salida VCO');
xlabel('Muestras');
ylabel('Amplitud');
