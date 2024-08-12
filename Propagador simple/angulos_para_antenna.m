% Devuelve el angulo entre el vector v1 y la proyeccion del vector v2 con
% el plano que definen v3 y v4
function [theta,flag] = angulos_para_antenna(v1,v2,v3,v4)

% Calcular el vector normal al plano definido por v3 y v4
normal = cross(v3, v4);
% Calcular la proyección de v2 sobre el plano
v2_proj = v2 - dot(v2, normal) / dot(normal,normal) * normal;
% a=-1*(dot(v2, normal) / dot(normal,normal) * normal);
% quiver3(0,0,0,a(1),a(2),a(3), 'LineWidth', 2,'Color','green') % eje x
% quiver3(0,0,0,v2_proj(1),v2_proj(2),v2_proj(3), 'LineWidth', 2,'Color','green') % eje x
% Calcular el ángulo entre v1 y v2_proj
cos_theta = dot(v1, v2_proj) / (norm(v1) * norm(v2_proj));
theta = acosd(cos_theta);

% Para salvar el error de ángulos unicamente positivos analizamos el ángulo
% entre el vector que queremos reflejar y el otro vector que corresponde a
% la definicion del plano (sabiendo que el plano está definido por la dirección de
% máxima radiación de la antena y este dicho vector, para nosotros v4) si
% este ángulo da mayor a 90 diremos que la componente ángular que estemos
% analizando es negativa, unicamente invertimos el signo por fuera de esta
% función, pero damos aviso poniendo una flag en 1.

ang_flag = acosd(dot(v2,v4)/(norm(v2)*norm(v4)));

if (ang_flag>90)
    flag = 1;
else
    flag = 0;
end

end