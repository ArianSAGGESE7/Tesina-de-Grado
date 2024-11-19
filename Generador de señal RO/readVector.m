function reconstructed_vector = readVector(output_folder)
    % Reconstruye el vector original leyendo archivos de pedazos guardados
    %
    % Inputs:
    %   output_folder: Carpeta donde están los archivos .mat con los pedazos
    %
    % Outputs:
    %   reconstructed_vector: Vector completo reconstruido a partir de los archivos

    % Obtener lista de archivos .mat en la carpeta
    file_list = dir(fullfile(output_folder, 'chunk_*.mat'));
    if isempty(file_list)
        error('No se encontraron archivos en la carpeta "%s".', output_folder);
    end

    % Leer los archivos y reconstruir el vector
    reconstructed_vector = [];
    for i = 1:length(file_list)
        % Cargar cada archivo
        data = load(fullfile(output_folder, file_list(i).name));
        reconstructed_vector = [reconstructed_vector, data.chunk];
    end

    fprintf('Vector reconstruido desde %d archivos.\n', length(file_list));
end
