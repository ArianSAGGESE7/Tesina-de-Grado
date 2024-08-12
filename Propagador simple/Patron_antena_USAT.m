%% Patrón de radiación de la antena

clc;clear;close all
Radiacion_antena1 = csvread('USAT-GNSS-L1-dB-RealizedGain.csv');
% Radiacion_antena2 = csvread('USAT-GNSS-L2-dB-RealizedGain.csv');
% Resolución de dos 6 grados
azimuth1 = -180:2:180;
elevation1 = 180:-2:0;
% patternCustom(Radiacion_antena1,azimuth1,elevation1,'Slice','theta','SliceValue',[0 180])
% patternCustom(Radiacion_antena1,azimuth1,elevation1)
% figure
% patternCustom(Radiacion_antena1,azimuth1,elevation1,CoordinateSystem="rectangular",Slice="theta",SliceValue=90)

% mesh(azimuth1,elevation1,Radiacion_antena1)

%% Reacomodamos filas y columnas de manera lineal
% phi azimut
% theta elevation
figure
tam=size(Radiacion_antena1);

Gain_antenna = reshape(Radiacion_antena1,[],1);

for i=1:tam(2)
    aux_matrix1(i,:) = azimuth1(i)*ones(1,tam(1));
end

azimuth = reshape(aux_matrix1',[],1);

for i=1:tam(2)
    aux_matrix(i,:) = elevation1;
end

elevation = reshape(aux_matrix',[],1);

% Patrón de antena completo
% patternCustom(Gain_antenna,elevation,azimuth)
 
A(:,1)=elevation;
A(:,2) = azimuth;
A(:,3)= Gain_antenna;

for i=1:length(Gain_antenna)

    if Gain_antenna(i)<3
        Gain_antenna_3dB(i) = 0;
    else
        Gain_antenna_3dB(i) = Gain_antenna(i);
    end
end
close all
% Patrón de antena de 3dB
figure(1)
patternCustom(Gain_antenna_3dB,elevation,azimuth)
figure(2)
patternCustom(Gain_antenna,elevation,azimuth)
% hacemos cortes para definir los ángulos de azimuth y elevacion
figure(3);
patternCustom(Gain_antenna_3dB,elevation,azimuth,CoordinateSystem="rectangular",Slice="phi",SliceValue=90)
figure(4);
patternCustom(Gain_antenna_3dB,elevation,azimuth,CoordinateSystem="rectangular",...
    Slice="phi ", SliceValue=0);
