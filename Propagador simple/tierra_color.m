function []=tierra_color()

space_color = 'w';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 0.3; % globe transparency level, 1 = opaque, through 0 = invisible
GMST0 = 0;%4.89496121282306; % Set up a rotatable globe at J2000.0
% Earth texture image
% image.
% image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
image_file = 'https://img.freepik.com/vector-gratis/vista-superior-fondo-mapa-mundo_1308-68322.jpg?w=1060&t=st=1709563009~exp=1709563609~hmac=9518fb2deb3ed318031153a257d1907f04082ba8af181d9316a036bdabaa33fd';
% Mean spherical earth
erad    = 6371008.7714; % equatorial radius (meters)
prad    = 6371008.7714; % polar radius (meters)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
% Create figure
figure('Color', space_color);
hold on;
% Turn off the normal axes
set(gca, 'NextPlot','add', 'Visible','off');
axis equal;
axis auto;
% Set initial view
view(0,30);
axis vis3d;
% Create wireframe globe
% Create a 3D meshgrid of the sphere points using the ellipsoid function
E = referenceEllipsoid(7030); % 7030 corresponde a WGS84
[x1, y1, z1] = ellipsoid(0,0,0,E.SemimajorAxis, E.SemimajorAxis, E.SemiminorAxis,180);
globe = surf(x1, y1, -z1, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end
% Texturemap the globe
% Load Earth image for texture map
cdata = imread(image_file);
% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha',alpha, 'EdgeColor', 'none');

% hold on
% theta = linspace(0, 2*pi, 100);
% x_equator = E.SemimajorAxis * cos(theta);
% y_equator = E.SemimajorAxis * sin(theta);
% z_equator = zeros(size(x_equator));
% plot3(x_equator, y_equator, z_equator, 'r', 'LineWidth', 2);


end