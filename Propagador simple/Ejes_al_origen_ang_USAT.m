% figure(2)
%% EJES EN EL ORIGEN

% Llevamos los ejes al origen 
pos=1;
% Eje desde el centro de la tierra al LEO 

a1 = -1*Posicion_LEO(:,pos);

a1 = a1/norm(a1); % Lo hacemos unitario 

% Eje perpendicular determinado por la velocidad 

a2 = Velocidad_LEO(:,pos);
a2 = a2/norm(a2);

% Eje perpendicular a la velocidad y la posición 

a3 = cross(a1,a2);
a3 = a3/norm(a3);

% Dibujamos los ejes
% tierra_color();
hold on

quiver3(0,0,0,a1(1)*10e6,a1(2)*10e6,a1(3)*10e6, 'LineWidth', 2,'Color','blue') % eje x
quiver3(0,0,0,a2(1)*10e6,a2(2)*10e6,a2(3)*10e6, 'LineWidth', 2,'Color','red') % eje y 
quiver3(0,0,0,a3(1)*10e6,a3(2)*10e6,a3(3)*10e6, 'LineWidth', 2,'Color','black') % eje z

quiver3(0,0,0,-1*y(1,end)*10e6 ,-1*y(2,end)*10e6 ,-1*y(3,end)*10e6 ,'LineWidth', 2,'Color','magenta')

%% CÁLCULO DE ÁNGULOS ENTRE LA ÚLTIMA DIRECCIÓN INICIAL Y LOS EJES QUE DEFINIMOS 

% Los ángulos en particular que queremos es entre -y0 y los ejes 
% Podemos calcular los ángulos sin necesidad de hacer una translación al
% origen, aunque podriamos verlo ahí. 

cos_dir_x = dot(a1,y(:,end))/(norm(a1));

cos_dir_z = dot(a3,y(:,end))/(norm(a3));

% Ahora analizamos si estos angulos caben dentro del +/- 50 grados respecto
% del eje "y" en la dirección "x" , y +/- 20 grados en el eje "y" en la dirección de "z" 

ang_x = acos(cos_dir_x)*180/pi-90;
ang_z = acos(cos_dir_z)*180/pi; 
