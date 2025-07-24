function plot_population_psth_from_multiple_structs(ob_structs, pcx_structs)
% Plot population PSTH with SEM for OB and PCx using multiple input structs

    edges   = -0.1:0.01:2;           % 10 ms bins from -100 to +400 ms
    t       = edges(1:end-1) + 0.005;  % bin centers
    binSize = diff(edges(1:2));        % seconds

    % ----- Merge OB and PCx safely -----
    ob_all  = cat_structs(ob_structs);
    pcx_all = cat_structs(pcx_structs);

    % ----- Compute PSTH matrices -----
    ob_matrix  = get_psth_matrix(ob_all, edges, binSize);
    pcx_matrix = get_psth_matrix(pcx_all, edges, binSize);

    % ----- Compute mean and SEM -----
    ob_mean = mean(ob_matrix, 1);
    ob_sem  = std(ob_matrix, 0, 1) / sqrt(size(ob_matrix, 1));

    pcx_mean = mean(pcx_matrix, 1);
    pcx_sem  = std(pcx_matrix, 0, 1) / sqrt(size(pcx_matrix, 1));

    % ----- Plot -----
    figure; hold on;

    % OB trace
    fill([t fliplr(t)], [ob_mean + ob_sem, fliplr(ob_mean - ob_sem)], ...
        [0.8 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(t, ob_mean, 'r', 'LineWidth', 2, 'DisplayName', 'OB');

    % PCx trace
    fill([t fliplr(t)], [pcx_mean + pcx_sem, fliplr(pcx_mean - pcx_sem)], ...
        [0.2 0.2 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(t, pcx_mean, 'k', 'LineWidth', 2, 'DisplayName', 'PCx');

    % Shaded response windows
    yMax = max([ob_mean + ob_sem, pcx_mean + pcx_sem]) * 1.1;
    patch([0 0.1 0.1 0], [0 0 yMax yMax], [0.9 0.9 0.9], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);
    patch([0.1 0.3 0.3 0.1], [0 0 yMax yMax], [0.7 0.7 0.7], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Formatting
    xline(0, '--k');
    xlabel('Time (s)');
    ylabel('Firing Rate (Hz)');
    title('Population PSTH: OB vs. PCx');
    legend('show');
    xlim([min(t), max(t)]);
    ylim([0, yMax]);
    set(gca, 'FontSize', 12);
end

% ===================== Robust Struct Concatenation ===================== %

function merged = cat_structs(structCell)
% Robustly concatenate structs with differing fields or field order
% Excludes 'vials' from tabular concatenation and reattaches afterward

    all_fields = {};
    has_vials = false;

    for i = 1:numel(structCell)
        if isempty(structCell{i}), continue; end
        all_fields = union(all_fields, fieldnames(structCell{i}));
    end

    % Remove 'vials' temporarily
    all_fields_no_vials = setdiff(all_fields, {'vials'});
    has_vials = any(strcmp(all_fields, 'vials'));

    all_tables = {};
    vials_all  = {};

    for i = 1:length(structCell)
        s = structCell{i};
        if isempty(s), continue; end

        nStructs = numel(s);  % s might be a struct array

        for j = 1:nStructs
            current = s(j);

            % Save and strip vials
            if has_vials && isfield(current, 'vials')
                vials_all{end+1,1} = current.vials;
                current = rmfield(current, 'vials');
            else
                vials_all{end+1,1} = [];
            end

            all_tables{end+1,1} = struct_to_padded_table(current, all_fields_no_vials);
        end
    end

    big_table = vertcat(all_tables{:});
    merged = table2struct(big_table);

    % Reattach vials
    if has_vials
        for i = 1:length(merged)
            merged(i).vials = vials_all{i};
        end
    end
end

function T = struct_to_padded_table(s, allFields)
% Pads struct s to contain all fields in allFields, then converts to table (row-wise)

    s_fields = fieldnames(s);
    missing = setdiff(allFields, s_fields);
    for i = 1:numel(missing)
        [s.(missing{i})] = deal([]);
    end
    s = orderfields(s, allFields);  % enforce same field order

    % Use AsArray = true to support scalar struct with fields of varying sizes
    T = struct2table(s, 'AsArray', true);
end

% ===================== PSTH Matrix Calculation ===================== %

function rate_matrix = get_psth_matrix(data, edges, binSize)
% Return matrix [units × bins] of avg firing rates per unit

    nUnits = numel(data);
    nBins = length(edges) - 1;
    rate_matrix = zeros(nUnits, nBins);

    for i = 1:nUnits
        spikes_all = [];

        for v = 1:numel(data(i).vials)
            trials = data(i).vials(v).trials;
            if isempty(trials), continue; end

            for t = 1:min(10, numel(trials))
                st = trials{t}(:);  % ensure column
                st = st(st >= edges(1) & st <= edges(end));
                spikes_all = [spikes_all; st]; %#ok<AGROW>
            end
        end

        if ~isempty(spikes_all)
            counts = histcounts(spikes_all, edges);
            rate_matrix(i,:) = counts / (binSize * 10); % Hz
        end
    end
end
