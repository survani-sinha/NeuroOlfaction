function process_cousens_units_and_append()
    % Load metadata
    meta = readtable('cousens.csv', 'VariableNamingRule', 'preserve');

    % Define binning for post-odor inhibition window
    post_window = [0 2];
    bin_size = 0.05;
    edges = post_window(1):bin_size:post_window(2);
    bin_centers = edges(1:end-1);

    inhibitionData = {};

    for i = 1:height(meta)
        % Clean variable name
        raw_name = regexprep(meta.("Session Name"){i}, '\s', '');  % removes all whitespace
        varname = regexprep(raw_name, '\..*', '');  % remove .00
        region = strtrim(meta.("Anatomical Location"){i});
        unit = meta.Unit(i);
        baseline_rate = meta.("Basal (pre-odor) Firing Rate (Hz)")(i);

        % Check existence in base workspace
        if ~evalin('base', sprintf("exist('%s', 'var')", varname))
            warning('⚠️ Variable "%s" not found in workspace.', varname);
            continue;
        end

        data = evalin('base', varname);

        % Validate .vials
        if ~isfield(data, 'vials') || ~iscell({data.vials})
            warning('⚠️ "%s" is missing valid .vials field or structure.', varname);
            continue;
        end

        vials = data.vials;
        if ~isstruct(vials)
            warning('⚠️ "%s.vials" is not a struct array.', varname);
            continue;
        end

        for vial = 1:numel(vials)
            if ~isfield(vials(vial), 'trials') || isempty(vials(vial).trials)
                continue;
            end

            trials = vials(vial).trials;

            for t = 1:size(trials, 2)
                spikes = trials{1, t};
                if isempty(spikes)
                    continue;
                end

                % Histogram firing rate over post-odor window
                counts = histcounts(spikes, edges);
                rates = counts / bin_size;

                % First time point where rate < baseline
                below_idx = find(rates < baseline_rate, 1);
                if ~isempty(below_idx)
                    inhibit_time = bin_centers(below_idx);
                    if inhibit_time < 0.05
                        continue;
                    end
                else
                    continue;
                end

                % Append result
                inhibitionData = [inhibitionData;
                    {varname, unit, vial, region, inhibit_time}];
            end
        end
    end

    % Convert to table
    T = cell2table(inhibitionData, 'VariableNames', ...
        {'File', 'Unit', 'Vial', 'Region', 'InhibitionOnsetTime'});

    % Append to existing file
    if isfile('inhibition_onset_summary.csv')
        existing = readtable('inhibition_onset_summary.csv');
        combined = [existing; T];
    else
        combined = T;
    end

    writetable(combined, 'inhibition_onset_summary.csv');
    disp('✅ Cousens data successfully added to inhibition_onset_summary.csv');
end
