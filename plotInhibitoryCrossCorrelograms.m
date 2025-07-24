function plotInhibitoryCrossCorrelograms(data)
% plotInhibitoryCrossCorrelograms(data)
%   Scans all cell pairs for putative inhibition (CCH troughs post-zero).
%   Plots only inhibitory pairs and writes summary CSV.

    numCells = numel(data);
    allPairs = nchoosek(1:numCells, 2);
    numTotalPairs = size(allPairs, 1);

    fprintf('Scanning %d unique cell pairs for inhibition...\n', numTotalPairs);

    % Initialize summary
    summary = [];
    sigCount = 0;

    for i = 1:numTotalPairs
        cellA = allPairs(i,1);
        cellB = allPairs(i,2);
        [isInhib, peakLag, minCount, baseline, baselineSD] = ...
            checkInhibitoryConnection(data, cellA, cellB, i, numTotalPairs);

        if isInhib
            sigCount = sigCount + 1;
            summary(sigCount).CellA = cellA;
            summary(sigCount).CellB = cellB;
            summary(sigCount).TroughLag = peakLag;
            summary(sigCount).MinCount = minCount;
            summary(sigCount).Baseline = baseline;
            summary(sigCount).BaselineSD = baselineSD;
        end
        close(gcf);
    end

    if sigCount > 0
        summaryTable = struct2table(summary);
        writetable(summaryTable, 'cross_correlogram_inhibition_summary.csv');
        fprintf('✓ Found %d inhibitory pairs. Summary saved.\n', sigCount);
    else
        fprintf('No inhibitory connections detected.\n');
    end
end

function [isInhibitory, troughLag, minCount, baseline, baselineSD] = ...
         checkInhibitoryConnection(data, cellA, cellB, idx, total)

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
    shiftCountsAll = zeros(numShuffles, length(centers));
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

    % Threshold
    nBins = length(centers);
    flankWidth = round(nBins / 4);
    baselineBins = [1:flankWidth, (nBins - flankWidth + 1):nBins];
    baseline = mean(realCounts(baselineBins));
    baselineSD = std(shiftCountsAll(:,baselineBins), 0, 'all');
    lowerThresh = baseline - SD_THRESHOLD * baselineSD;

    % Find inhibitory bins *after* time = 0
    postZero = centers > 0;
    sigLow = realCounts < lowerThresh;
    inhibBins = sigLow & postZero;

    % Define as inhibitory if >= 3 bins below threshold after 0
    [minCount, minIdx] = min(realCounts);
    troughLag = centers(minIdx);
    isInhibitory = sum(inhibBins) >= 3;

    if isInhibitory
        % Plot
        figure('Name', sprintf('INHIBITION: Cell %d → %d', cellA, cellB));
        hold on;
        bar(centers, realCounts, 'k');
        yline(baseline, 'g', 'LineWidth', 1);
        yline(lowerThresh, 'b--', 'LineWidth', 1);
        scatter(centers(inhibBins), realCounts(inhibBins), 25, 'c', 'filled');
        xlabel('Time lag (s)'); ylabel('Spike count');
        title(sprintf('INHIBITORY CCH: Cell %d → %d  (%d/%d)', cellA, cellB, idx, total));
        legend('Real', 'Baseline', 'Lower Thresh', 'Inhibitory Bins');
        xlim([-maxLag maxLag]);
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
