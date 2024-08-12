%% Transformación al sistema de coordenadas en el USAT 

function [T_org2usat] = transformacion_USAT(pos_USAT,vel_USAT)
% Primero definimos los vectores que dan lugar a nuestro nuevo sistema de
% referencia 

axis_x_usat = -pos_USAT/norm(pos_USAT);

axis_y_usat = vel_USAT/norm(vel_USAT);

axis_z_usat = cross(axis_y_usat,axis_x_usat);

axis_z_usat = axis_z_usat/(norm(axis_z_usat));

T_org2usat = [axis_x_usat axis_y_usat axis_z_usat];

end
