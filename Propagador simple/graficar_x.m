function []=graficar_x(x)

[filas,columnas] =size(x);

for i = 1:filas

    plot3([x(i,1)],[x(i,2)],[x(i,3)], 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'black','LineWidth',2)

end

end