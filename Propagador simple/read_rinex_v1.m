function [header, observations, keplerian_params] = read_rinex_v1(filename)
    fid = fopen(filename, 'r');
    if fid == -1
        error('Unable to open file');
    end
    
    % Read header
    header = cell(0);
    while true
        line = fgetl(fid);
        if ~ischar(line) || strncmp(line, 'END OF HEADER', 13)
            break;
        end
        header{end+1} = line;
    end
    
    % Read observation data
    observations = cell(0);
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        observations{end+1} = line;
    end
    
    % Read navigation data and extract Keplerian parameters
    keplerian_params = struct('prn', [], 'semimajor_axis', [], 'eccentricity', [], ...
                              'inclination', [], 'argument_perigee', [], ...
                              'longitude_asc_node', [], 'mean_anomaly', []);
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        disp(line); % Display the line to check its content
        if strncmp(line, 'NAV', 3) % Check if it's a navigation data line
            prn = str2double(line(4:6)); % PRN number
            semimajor_axis = str2double(line(24:42)); % Semimajor axis (m)
            eccentricity = str2double(line(43:61)); % Eccentricity
            inclination = str2double(line(62:80)); % Inclination (degrees)
            argument_perigee = str2double(line(81:99)); % Argument of perigee (degrees)
            longitude_asc_node = str2double(line(100:118)); % Longitude of ascending node (degrees)
            mean_anomaly = str2double(line(119:137)); % Mean anomaly (degrees)
            
            % Store Keplerian parameters in struct array
            keplerian_params(end+1) = struct('prn', prn, 'semimajor_axis', semimajor_axis, ...
                                              'eccentricity', eccentricity, 'inclination', inclination, ...
                                              'argument_perigee', argument_perigee, ...
                                              'longitude_asc_node', longitude_asc_node, ...
                                              'mean_anomaly', mean_anomaly);
        end
    end
    
    % Remove the empty struct at the beginning
    keplerian_params = keplerian_params(2:end);
    
    fclose(fid);
end
