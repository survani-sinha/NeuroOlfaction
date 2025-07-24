function browseAutocorrelogramsByCell(dataInput)
% browseAutocorrelogramsByCell(dataInput)
%   Computes and browses autocorrelograms for each unit in the data.
%   Groups by unit index (cell), assuming each data(u) is a cell.
%   Use ← and → keys to navigate between cells.
%
%   Example usage:
%       browseAutocorrelogramsByCell('NEW_aPCx_*')
%       browseAutocorrelogramsByCell(NEW_aPCx_15)

    binSize = 0.001;     % 1 ms
    maxLag = 0.5;        % +/- 500 ms
    lags = -maxLag:binSize:maxLag;

    allCellACGs = {};    % stores autocorr results per cell

    % Load data
    if ischar(dataInput) || isstring(dataInput)
        varNames = who(dataInput);
        if isempty(varNames)
            error('No variables match pattern "%s".', dataInput);
        end
        dataList = cellfun(@(n) evalin('base', n), varNames, 'UniformOutput', false);
    elseif isstruct(dataInput)
        dataList = {dataInput};
    else
        error('Input must be string or struct.');
    end

    % For each cell (data(u)), gather all spike times
    for d = 1:length(dataList)
        data = dataList{d};
        for u = 1:numel(data)
            allSpikes = [];
            for vial = 1:numel(data(u).vials)
                trials = data(u).vials(vial).trials;
                for t = 1:numel(trials)
                    st = trials{t};
                    allSpikes = [allSpikes; st(:)];
                end
            end
            allSpikes = sort(allSpikes);
            if numel(allSpikes) > 100
                acg = computeAutocorr(allSpikes, binSize, maxLag);
                allCellACGs{end+1} = acg;
            else
                allCellACGs{end+1} = zeros(size(lags));
            end
        end
    end

    currentCell = 1;
    hFig = figure('Name','Autocorrelogram Browser','KeyPressFcn',@onKey);
    drawCell();

    function onKey(~, evt)
        switch evt.Key
            case 'rightarrow'
                currentCell = min(currentCell + 1, numel(allCellACGs));
            case 'leftarrow'
                currentCell = max(currentCell - 1, 1);
            otherwise
                return;
        end
        drawCell();
    end

    function drawCell()
        clf(hFig);
        bar(lags, allCellACGs{currentCell}, 'k', 'BarWidth', 1);
        xlim([-maxLag maxLag]);
        xlabel('Lag (s)');
        ylabel('Autocorrelation');
        title(sprintf('Autocorrelogram for Cell %d/%d', currentCell, numel(allCellACGs)));
    end
end

function acg = computeAutocorr(spikeTimes, binSize, maxLag)
    lags = -maxLag:binSize:maxLag;
    acg = zeros(size(lags));
    for i = 1:length(spikeTimes)
        diffs = spikeTimes - spikeTimes(i);
        diffs(i) = [];  % remove zero lag (self-count)
        diffs = diffs(abs(diffs) <= maxLag);
        bins = round((diffs + maxLag) / binSize) + 1;
        for b = bins(:)'
            acg(b) = acg(b) + 1;
        end
    end
end
