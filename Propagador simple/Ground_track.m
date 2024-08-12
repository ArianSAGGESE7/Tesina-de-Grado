function [] = Ground_track(sat_pos)
% Función para hacer el trackeo del satélite en la tierra, debe estar dado
% en ECEF

a=size(sat_pos);
for i=1:a(2)
    [lla1] = ecef2llaGeod(sat_pos(:,i)); 
    lla(i,:) = lla1;
end


% geoplot(lla(:,1)',lla(:,2)','-*')
hold on
worldmap('World')
load coastlines
plotm(coastlat,coastlon)

plotm(lla(:,1)',lla(:,2)','o','Color','red','Marker','+','LineWidth',2)


end


%% De 0 - 360 a -180 - 180

% lon0 = 200;
% lon1 = mod((lon0 - 180), 360) - 180
% 
% %% De -180 - 180 a 0 - 360
% 
% lon0 = -20;
% lon1 = mod(lon0, 360)