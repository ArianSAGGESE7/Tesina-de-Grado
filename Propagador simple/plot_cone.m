
function [x_circle,y_circle,z_circle] = plot_cone(P, V, R)

    % P: 1x3 array representing the center point of the circle
    % V: 1x3 array representing the vector defining the normal plane
    % R: radius of the circle
    
    % Normalize the vector V
    V = V / norm(V);
    
    % Define an arbitrary vector perpendicular to V
    if V(1) == 0 && V(2) == 0
        perp_vector = [1, 0, 0]; % Choose an arbitrary perpendicular vector
    else
        perp_vector = [V(2), -V(1), 0];
    end
    
    % Normalize the perpendicular vector
    perp_vector = perp_vector / norm(perp_vector);
    
    % Calculate the third orthogonal vector
    third_orthogonal_vector = cross(V, perp_vector);
    
    % Generate points along the circumference of the circle
    t = linspace(0, 2*pi, 100);
    
    % Calculate points on the circle
    x_circle = P(1) + R * (perp_vector(1) * cos(t) + third_orthogonal_vector(1) * sin(t));
    y_circle = P(2) + R * (perp_vector(2) * cos(t) + third_orthogonal_vector(2) * sin(t));
    z_circle = P(3) + R * (perp_vector(3) * cos(t) + third_orthogonal_vector(3) * sin(t));
end
