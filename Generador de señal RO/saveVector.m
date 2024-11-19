function saveVector(long_vector, chunk_size, output_folder)
    % Divide un vector largo en pedazos y los guarda en archivos .mat
    %
    % Inputs:
    %   long_vector: Vector largo que se quiere dividir
    %   chunk_size: Número máximo de elementos en cada archivo
    %   output_folder: Carpeta donde se guardarán los archivos

    % Crear la carpeta de salida si no existe
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % Dividir el vector en pedazos y guardar
    num_chunks = ceil(length(long_vector) / chunk_size); % Número total de pedazos
    for i = 1:num_chunks
        % Calcular los índices para el pedazo actual
        start_idx = (i - 1) * chunk_size + 1;
        end_idx = min(i * chunk_size, length(long_vector));
        
        % Extraer el pedazo actual
        chunk = long_vector(start_idx:end_idx);
        
        % Guardar el pedazo en un archivo .mat
        filename = fullfile(output_folder, sprintf('chunk_%d.mat', i));
        save(filename, 'chunk');
    end

    fprintf('Vector dividido y guardado en %d archivos en la carpeta "%s".\n', num_chunks, output_folder);
end
