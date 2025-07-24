function statsTable = computePSTHWilcoxonStats(dataInput)
% computePSTHWilcoxonStats(dataInput)
%   Performs Wilcoxon signed-rank test comparing pre-odor and during-odor
%   firing rates for each unit and each vial.
%
%   Input:
%     - A string pattern (e.g. 'NEW_aPCx_*') or a struct variable
%
%   Output:
%     - statsTable: table with file, unit, vial, stats & direction
%
%   Odor window = [0 2], Pre window = [-2 0]

    PST = [-4 8]; binWidth = 0.1;
    edges = PST(1):binWidth:PST(2);
    preWindow   = [-2 0];
    odorWindow  = [0 2];

    preIdx  = edges(1:end-1) >= preWindow(1)  & edges(1:end-1) < preWindow(2);
    odorIdx = edges(1:end-1) >= odorWindow(1) & edges(1:end-1) < odorWindow(2);

    stats = {};

    % Load data
    if ischar(dataInput) || isstring(dataInput)
        inputName = char(dataInput);
        if evalin('base', sprintf('exist(''%s'', ''var'')', inputName))
            dataList = {evalin('base', inputName)};
            fileList = {inputName};
        else
            varNames = who(inputName);
            if isempty(varNames)
                error('No variables match pattern or name "%s".', inputName);
            end
            dataList = cellfun(@(n) evalin('base', n), varNames, 'UniformOutput', false);
            fileList = varNames;
        end
    elseif isstruct(dataInput)
        dataList = {dataInput};
        fileList = {'(input struct)'};
    else
        error('Input must be a string or struct.');
    end

    for d = 1:numel(dataList)
        data = dataList{d};
        fname = fileList{d};

        for u = 1:numel(data)
            for vial = 1:numel(data(u).vials)
                trials = data(u).vials(vial).trials;
                allRates = [];

                for t = 1:numel(trials)
                    st = trials{t}(:);
                    st = st(st >= PST(1) & st <= PST(2));
                    rates = histcounts(st, edges) / binWidth;

                    if numel(rates) == numel(edges)-1
                        allRates(end+1,:) = rates; %#ok<AGROW>
                    end
                end

                if isempty(allRates)
                    continue;
                end

                preRates  = mean(allRates(:,preIdx), 2);
                odorRates = mean(allRates(:,odorIdx), 2);

                try
                    p = signrank(preRates, odorRates);
                catch
                    p = NaN;
                end

                if mean(odorRates) < mean(preRates)
                    direction = 'Inhibitory';
                elseif mean(odorRates) > mean(preRates)
                    direction = 'Excitatory';
                else
                    direction = 'None';
                end

                stats(end+1, :) = {
                    fname, u, vial, ...
                    mean(preRates), std(preRates), ...
                    mean(odorRates), std(odorRates), ...
                    p, direction
                }; %#ok<AGROW>
            end
        end
    end

    statsTable = cell2table(stats, 'VariableNames', {
        'File', 'Unit', 'Vial', ...
        'PreRateMean', 'PreRateStd', ...
        'OdorRateMean', 'OdorRateStd', ...
        'pValue', 'Direction'
    });

    disp(statsTable);
end