function inhibitionTable = computeInhibitionOnsetTimes()
% computeInhibitionOnsetTimes()
% Combines raw data + CSV annotations to estimate inhibition onset.
% Returns a table with:
% - File, Unit, Vial, Region, InhibitionOnsetTime (s)

    % Settings
    PST = [-4 8];
    binWidth = 0.1;
    edges = PST(1):binWidth:PST(2);
    centers = edges(1:end-1) + binWidth/2;
    preWindow = [-2 0];
    postWindow = [0 2];
    preIdx = edges(1:end-1) >= preWindow(1) & edges(1:end-1) < preWindow(2);
    postIdx = edges(1:end-1) >= postWindow(1) & edges(1:end-1) < postWindow(2);

    % Data file list
    dataFiles = {
        'NEW_aPCx_15', 'NEW_plCoA_15', ...
        'NEW_141209_bank1', 'NEW_141209_bank2', 'NEW_160819_bank1', 'NEW_160819_bank2', ...
        'NEW_160820_bank1', 'NEW_160820_bank2', 'NEW_170608_bank1', 'NEW_170608_bank2', ...
        'AON1110451', 'AON1120395', 'AON1120433', 'AON2310381', 'AON2311111', ...
        'AON310482', 'AON3314460', 'AON3711000', 'AON3711030', 'AON3711100', ...
        'AON410345', 'AON410384', 'AON440315', 'AON440350', 'AON440369', ...
        'DTT1010345', 'DTT110302', 'DTT210300', 'DTT210463', 'DTT210485', ...
        'DTT210510', 'DTT2611111', 'DTT2710365', 'DTT3214100', 'DTT3214190', ...
        'DTT3423609', 'DTT3812510', 'DTT3812535', 'DTT4012471', 'DTT4012475', ...
        'DTT4312407', 'DTT510472', 'DTT520557', 'DTT530328', 'DTT540295', ...
        'DTT731111', 'DTT910353', 'DTT910363', 'DTT910388', ...
        'VTT2612222', 'VTT2613333', 'VTT2710381', 'VTT2710397', 'VTT2710415', ...
        'VTT2710435', 'VTT3114160', 'VTT3124165', 'VTT610396', 'VTT610427', 'VTT620387'
    };

    inhibitionData = {};  % will hold rows: {file, unit, vial, region, onset time}

    for i = 1:numel(dataFiles)
        fname = dataFiles{i};
        if ~evalin('base', sprintf("exist('%s', 'var')", fname))
            warning('⚠️ Variable "%s" not found in base workspace.', fname);
            continue;
        end

        data = evalin('base', fname);
        cleanName = erase(fname, 'NEW_');  % Match CSV name formatting
        csvName = sprintf('FilteredWilcoxon - %s.csv', cleanName);

        if ~isfile(csvName)
            warning('⚠️ CSV "%s" not found.', csvName);
            continue;
        end

        csvInfo = readtable(csvName, 'VariableNamingRule', 'preserve');

        % Identify anatomy/region column
        varnames = csvInfo.Properties.VariableNames;
        anatomyCol = '';
        for j = 1:numel(varnames)
            nameLower = lower(strrep(varnames{j}, '_', ' '));
            if contains(nameLower, 'anatomical') || ...
               contains(nameLower, 'region') || ...
               contains(nameLower, 'location') || ...
               contains(nameLower, 'cortical')
                anatomyCol = varnames{j};
                break;
            end
        end
        if isempty(anatomyCol)
            warning('⚠️ No anatomical/cortical/location/region column found in "%s".', csvName);
            continue;
        end

        % Process each unit × vial
        for u = 1:numel(data)
            for vial = 1:numel(data(u).vials)
                trials = data(u).vials(vial).trials;
                if isempty(trials)
                    continue;
                end

                try
                    matchRows = strcmp(strtrim(string(csvInfo.File)), [fname '''']) & ...
                                csvInfo.Unit == u & ...
                                csvInfo.Vial == vial;
                catch
                    warning('⚠️ Could not match file "%s" in CSV.', fname);
                    continue;
                end

                if ~any(matchRows)
                    continue;
                end

                region = csvInfo.(anatomyCol){find(matchRows, 1)};
                if isempty(region) || strcmpi(strtrim(region), 'Unknown')
                    continue;
                end

                % PSTH calculation
                trialRates = zeros(numel(trials), numel(centers));
                for t = 1:numel(trials)
                    st = trials{t}(:);
                    st = st(st >= PST(1) & st <= PST(2));
                    trialRates(t,:) = histcounts(st, edges) / binWidth;
                end

                psth = mean(trialRates, 1);
                baseline = mean(psth(preIdx));
                postBins = find(postIdx);
                belowIdx = find(psth(postBins) < baseline, 1);

                if isempty(belowIdx)
                    continue;
                end

                inhibitTime = centers(postBins(belowIdx));
                inhibitionData(end+1, :) = {fname, u, vial, region, inhibitTime};
            end
        end
    end

    % Output table
    if isempty(inhibitionData)
        warning('⚠️ No inhibition data extracted.');
        inhibitionTable = table();
    else
        inhibitionTable = cell2table(inhibitionData, 'VariableNames', ...
            {'File', 'Unit', 'Vial', 'Region', 'InhibitionOnsetTime'});
        writetable(inhibitionTable, 'inhibition_onset_summary.csv');
        disp('✅ Saved to inhibition_onset_summary.csv (excluding Unknown regions and NaNs)');
    end
end
