function plotCrossCorrelogram(data)
% plotAllCrossCorrelograms(data)
%   Interactively browse significant cross-correlograms (excitation & inhibition),
%   output figures, and generate a summary CSV file.

    numCells = numel(data);
    allPairs = nchoosek(1:numCells, 2);
    numTotalPairs = size(allPairs, 1);

    fprintf('Found %d cells (%d unique pairs). Scanning for significance...\n', ...
            numCells, numTotalPairs);

    % Initialize summary table
    summary = [];

    useCached = false;

    if useCached && isfile('significant_pairs.mat')
        load('significant_pairs.mat', 'sigPairs');
        sigCount = size(sigPairs, 1);
        fprintf('→ Loaded %d significant pairs from cache.\n', sigCount);
    else
        sigPairs = [];
        sigCount = 0;
        for i = 1:numTotalPairs
            cellA = allPairs(i,1);
            cellB = allPairs(i,2);
            [isSig, peakLag, peakType, maxCount, baseline, baselineSD] = ...
                plotCrossCorrelogramIn(data, cellA, cellB, i, numTotalPairs, false);

            if isSig
                sigCount = sigCount + 1;
                sigPairs(sigCount,:) = [cellA, cellB]; %#ok<AGROW>

                % Append to summary
                summary(sigCount).CellA = cellA;
                summary(sigCount).CellB = cellB;
                summary(sigCount).PeakLag = peakLag;
                summary(sigCount).Type = peakType;
                summary(sigCount).MaxCount = maxCount;
                summary(sigCount).Baseline = baseline;
                summary(sigCount).BaselineSD = baselineSD;
            end
            close(gcf);
        end
        fprintf('→ Found %d significant pairs. Entering interactive mode...\n', sigCount);
        save('significant_pairs.mat', 'sigPairs', 'numCells', 'numTotalPairs');
    end

    % Save summary
    if sigCount > 0 && ~isempty(summary)
        summaryTable = struct2table(summary);
        writetable(summaryTable, 'cross_correlogram_summary.csv');
        fprintf('→ CSV summary written to cross_correlogram_summary.csv\n');
    else
        fprintf('→ No significant pairs to summarize.\n');
        return;
    end

    % INTERACTIVE BROWSER
    pairIdx = 1;
    hFig = figure('Name', 'Significant CCH Browser', 'NumberTitle', 'off', ...
                  'KeyPressFcn', @onKey);
    drawPlot();

    function onKey(~, evt)
        switch evt.Key
            case 'leftarrow'
                pairIdx = mod(pairIdx - 2, sigCount) + 1;
            case 'rightarrow'
                pairIdx = mod(pairIdx, sigCount) + 1;
            otherwise
                return;
        end
        drawPlot();
    end

    function drawPlot()
        clf(hFig);
        cellA = sigPairs(pairIdx, 1);
        cellB = sigPairs(pairIdx, 2);
        fprintf('→ Plotting SIGNIFICANT pair: Cell %d vs %d (%d/%d)\n', ...
                cellA, cellB, pairIdx, sigCount);
        plotCrossCorrelogramIn(data, cellA, cellB, pairIdx, sigCount, true);
    end
end

function [isSignificant, peakLag, peakType, maxCount, baseline, baselineSD] = ...
         plotCrossCorrelogramIn(data, cellA, cellB, idx, total, doPlot)

    if nargin < 6, doPlot = true; end

    maxLag = 0.1; binSize = 0.001;
    edges = -maxLag:binSize:maxLag;
    centers = edges(1:end-1) + binSize/2;
    numShuffles = 10;
    SD_THRESHOLD = 4.41;

    spikesA = getSpikeTimesByTrial(data(cellA));
    spikesB = getSpikeTimesByTrial(data(cellB));
    numTrials = min(numel(spikesA), numel(spikesB));

    % Real CCH
    realDiffs = [];
    for t = 1:numTrials
        a = spikesA{t}; b = spikesB{t};
        for i = 1:numel(a)
            dt = b - a(i);
            realDiffs = [realDiffs; dt(dt >= -maxLag & dt <= maxLag)];
        end
    end
    realCounts = histcounts(realDiffs, edges);

    % Shift predictor
    shiftCountsAll = zeros(numShuffles, length(edges)-1);
    for s = 1:numShuffles
        shuffledB = spikesB(randperm(numTrials));
        diffs = [];
        for t = 1:numTrials
            a = spikesA{t}; b = shuffledB{t};
            for i = 1:numel(a)
                dt = b - a(i);
                diffs = [diffs; dt(dt >= -maxLag & dt <= maxLag)];
            end
        end
        shiftCountsAll(s,:) = histcounts(diffs, edges);
    end

    % Significance
    flankWidth = round(length(centers) / 4);
    nBins = length(centers);  % ✅ FIXED: replaced use of `end` outside indexing
    baselineBins = [1:flankWidth, (nBins - flankWidth + 1):nBins];
    baseline = mean(realCounts(baselineBins));
    baselineSD = std(shiftCountsAll(:,baselineBins), 0, 'all');

    upperThresh = baseline + SD_THRESHOLD * baselineSD;
    lowerThresh = baseline - SD_THRESHOLD * baselineSD;

    sigHigh = realCounts > upperThresh;
    sigLow = realCounts < lowerThresh;

    [maxCount, maxIdx] = max(realCounts);
    peakLag = centers(maxIdx);

    isExcitatory = (sum(sigHigh) >= 3) && (maxCount > baseline + 6 * baselineSD);
    isInhibitory = sum(sigLow) >= 3;
    isSignificant = isExcitatory || isInhibitory;

    if isExcitatory && isInhibitory
        peakType = 'Both';
    elseif isExcitatory
        peakType = 'Excitatory';
    elseif isInhibitory
        peakType = 'Inhibitory';
    else
        peakType = 'None';
    end

    if doPlot
        hold on;
        bar(centers, realCounts, 'k');
        yline(baseline, 'g', 'LineWidth', 1);
        yline(upperThresh, 'r--', 'LineWidth', 1);
        yline(lowerThresh, 'b--', 'LineWidth', 1);
        scatter(centers(sigHigh), realCounts(sigHigh), 25, 'r', 'filled');
        scatter(centers(sigLow), realCounts(sigLow), 25, 'c', 'filled');
        xlim([-maxLag maxLag]);
        xlabel('Time lag (s)');
        ylabel('Spike count');
        title(sprintf('CCH: Cell %d vs %d (%d/%d)  %s', ...
            cellA, cellB, idx, total, ternary(isSignificant, '[Significant]', '')));
        legend('Real', 'Baseline', 'Upper thresh', 'Lower thresh', ...
               'Excitatory bins', 'Inhibitory bins');
        grid on;
    end
end

function spikeTrials = getSpikeTimesByTrial(unit)
    spikeTrials = {};
    for vial = 1:numel(unit.vials)
        trials = unit.vials(vial).trials;
        for t = 1:numel(trials)
            spikeTrials{end+1} = trials{t}(:); %#ok<AGROW>
        end
    end
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
