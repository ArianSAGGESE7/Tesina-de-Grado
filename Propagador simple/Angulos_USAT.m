%% Ángulo de llegada al USAT 

% Cálculo de la matriz de transformación ejes USAT 
[T_org2usat] = transformacion_USAT(Posicion_LEO(:,end),Velocidad_LEO(:,end));
hold on
quiver3(0, 0, 0,T_org2usat(1,1) , T_org2usat(2,1),T_org2usat(3,1) , 'LineWidth', 2,'Color','blue');
quiver3(0, 0, 0,T_org2usat(1,2) , T_org2usat(2,2),T_org2usat(3,2) , 'LineWidth', 2,'Color','red');
quiver3(0, 0, 0,T_org2usat(1,3) , T_org2usat(2,3),T_org2usat(3,3) , 'LineWidth', 2,'Color','black');
