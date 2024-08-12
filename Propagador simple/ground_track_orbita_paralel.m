% Script de simulación en paralelo de ground track + orbita
% Despues de simular Prueba_de_propagador podemos correr este script 

close all;
tierra_color() % simulamos el satélite y su órbita en este plot figura(1)
view(3)
figure(2) % figure(2) es el ground track
worldmap('World')
load coastlines
plotm(coastlat,coastlon)

%%
pause(10)
for i=1:40:length(pos_pos(:,1))
    figure(1)
    hold on
    plot3(pos_pos(i, 1), pos_pos(i, 2), pos_pos(i, 3),'Marker','diamond','Color','red')
    figure(2)
    hold on
    plotm(lla(i,1)',lla(i,2)'+90,'-o','Color','k')
end

