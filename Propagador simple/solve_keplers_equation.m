
function E = solve_keplers_equation(M, e)
    % Solución a la ecuación de Kepler usando el método de Newton 
    tolerance = 1e-12;
    E0 = M;
    diff = tolerance + 1;
    while diff > tolerance
        E = E0 - (E0 - e * sin(E0) - M) / (1 - e * cos(E0));
        diff = abs(E - E0);
        E0 = E;
    end
end
