% === CONFIGURATIONN ===
base_dir = 'C:\Users\Manuel Zúniga\Documents\Laboratorio\Imágenes\SEM\Análisis tetillas';
genotipos = {'Control', 'kay1', 'kay5', 'kayDF', 'Retinin'};
resumen = {};  % Para registrar los primeros picos

% === LOOP BY GENOTYPES Y EXPERIMENTS ===
for g = 1:length(genotipos)
    genotipo = genotipos{g};
    genotipo_path = fullfile(base_dir, genotipo);
    if ~isfolder(genotipo_path), continue; end

    experimento_folders = dir(genotipo_path);
    for e = 1:length(experimento_folders)
        carpeta = experimento_folders(e);
        if ~carpeta.isdir || startsWith(carpeta.name, '.'), continue; end

        carpeta_path = fullfile(genotipo_path, carpeta.name);
        archivo_csv = fullfile(carpeta_path, 'coordenadas.csv');

        if ~isfile(archivo_csv)
            fprintf('❌ No se encontró coordenadas.csv en %s\n', carpeta_path);
            continue;
        end

        try
            % === Read Data ===
            data = readtable(archivo_csv);
            X = data.x1;
            Y = data.y1;

            % ===  RDF Parameters===
            dr = 0.035;
            r_max = 1;
            d_min = 0.210;
            xmin = min(X); xmax = max(X);
            ymin = min(Y); ymax = max(Y);
            area = (xmax - xmin) * (ymax - ymin);
            rho = length(X) / area;

            % === Corrección por borde ===
            margin = r_max;
            valid_idx = (X > xmin + margin) & (X < xmax - margin) & ...
                        (Y > ymin + margin) & (Y < ymax - margin);
            X_valid = X(valid_idx); Y_valid = Y(valid_idx);
            N_valid = length(X_valid);

            edges = 0:dr:r_max;
            r = edges(1:end-1) + dr/2;
            g_r = zeros(size(r));

            for i = 1:N_valid
                dx = X - X_valid(i);
                dy = Y - Y_valid(i);
                d = sqrt(dx.^2 + dy.^2);
                d = d(d >= d_min);
                counts = histcounts(d, edges);
                shell_areas = pi * (edges(2:end).^2 - edges(1:end-1).^2);
                g_r = g_r + counts ./ shell_areas;
            end
            g_r = g_r / (rho * N_valid);

            % === Simulation ===
            N_rand = length(X);
            X_rand_hc = []; Y_rand_hc = [];
            attempts = 0; max_attempts = 1e5;
            while length(X_rand_hc) < N_rand && attempts < max_attempts
                x_try = xmin + (xmax - xmin) * rand;
                y_try = ymin + (ymax - ymin) * rand;
                if isempty(X_rand_hc)
                    accept = true;
                else
                    dists = sqrt((X_rand_hc - x_try).^2 + (Y_rand_hc - y_try).^2);
                    accept = all(dists >= d_min);
                end
                if accept
                    X_rand_hc(end+1) = x_try;
                    Y_rand_hc(end+1) = y_try;
                end
                attempts = attempts + 1;
            end
            X_rand = X_rand_hc(:); Y_rand = Y_rand_hc(:);
            valid_idx_rand = (X_rand > xmin + margin) & (X_rand < xmax - margin) & ...
                             (Y_rand > ymin + margin) & (Y_rand < ymax - margin);
            X_valid_rand = X_rand(valid_idx_rand);
            Y_valid_rand = Y_rand(valid_idx_rand);
            N_valid_rand = length(X_valid_rand);
            g_rand = zeros(size(r));
            for i = 1:N_valid_rand
                dx = X_rand - X_valid_rand(i);
                dy = Y_rand - Y_valid_rand(i);
                d = sqrt(dx.^2 + dy.^2);
                d = d(d >= d_min);
                counts = histcounts(d, edges);
                shell_areas = pi * (edges(2:end).^2 - edges(1:end-1).^2);
                g_rand = g_rand + counts ./ shell_areas;
            end
            g_rand = g_rand / (rho * N_valid_rand);

            % === First peak ===
            [peaks, locs] = findpeaks(g_r, r, 'MinPeakHeight', 1.05);
            if ~isempty(peaks)
                first_peak_pos = locs(1);
                first_peak_val = peaks(1);
            else
                first_peak_pos = NaN;
                first_peak_val = NaN;
            end

            % === Save data ===
            resultado_csv = table(r', g_r', g_rand', ...
                'VariableNames', {'r_um', 'g_r', 'g_r_rand'});
            writetable(resultado_csv, fullfile(carpeta_path, 'rdf_resultados.csv'));

            % === Save graph ===
            f = figure('Visible', 'off');
            plot(r, g_r, 'b-', 'LineWidth', 2); hold on;
            plot(r, g_rand, 'r--', 'LineWidth', 2);
            if ~isnan(first_peak_pos)
                plot(first_peak_pos, first_peak_val, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
                text(first_peak_pos, first_peak_val + 0.05, ...
                    sprintf('%.2f µm', first_peak_pos), ...
                    'HorizontalAlignment', 'center');
            end
            xlabel('Distancia r (\mum)');
            ylabel('g(r)');
            legend('Datos reales', 'Aleatoria (con exclusión)');
            title(sprintf('RDF - %s/%s', genotipo, carpeta.name));
            grid on;
            saveas(f, fullfile(carpeta_path, 'rdf_grafica.png'));
            close(f);

           
            resumen(end+1, :) = {genotipo, carpeta.name, first_peak_pos, first_peak_val};

        catch err
            fprintf('⚠️ Error procesando %s: %s\n', carpeta_path, err.message);
        end
    end
end

% === Save final table ===
resumen_tbl = cell2table(resumen, ...
    'VariableNames', {'Genotipo', 'Experimento', 'Pico_um', 'g_r_valor'});
writetable(resumen_tbl, fullfile(base_dir, 'resumen_RDF.csv'));