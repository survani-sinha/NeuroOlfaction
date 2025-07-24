function odorFiringHeatmap(data)
    odorPeriod = [0, 2];  % Odor window (in sec)
    nCells = numel(data);
    maxVials = max(cellfun(@(x) numel(x.vials), num2cell(data)));

    rateMatrix = NaN(nCells, maxVials);  % rows: cells, cols: vials

    for c = 1:nCells
        for v = 1:numel(data(c).vials)
            trials = data(c).vials(v).trials;
            rates = zeros(numel(trials), 1);

            for t = 1:numel(trials)
                spikes = trials{t};
                if isempty(spikes), continue; end
                rates(t) = sum(spikes >= odorPeriod(1) & spikes <= odorPeriod(2)) / diff(odorPeriod);
            end

            % average across trials for this vial
            rateMatrix(c, v) = mean(rates, 'omitnan');
        end
    end

    % Z-score normalize across vials per cell (row-wise)
    normMatrix = (rateMatrix - mean(rateMatrix, 2, 'omitnan')) ./ std(rateMatrix, 0, 2, 'omitnan');

    % Plot heatmap
    figure;
    imagesc(normMatrix);
    colormap('hot');
    colorbar;
    xlabel('Vial (Odorant #)');
    ylabel('Cell');
    title('Normalized Odor-Period Firing Rates');
end
