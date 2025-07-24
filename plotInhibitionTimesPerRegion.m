function plotInhibitionSummaries()
% plotInhibitionSummaries()
% Generates bar plot (with error bars) and boxplot of inhibition onset times by region.
% Requires 'inhibition_onset_summary.csv'

    %% --- Load Data ---
    summaryFile = 'inhibition_onset_summary.csv';
    if ~isfile(summaryFile)
        error('File "%s" not found.', summaryFile);
    end

    tbl = readtable(summaryFile);
    tbl = tbl(~isnan(tbl.InhibitionOnsetTime), :);  % remove NaNs

    regions = unique(tbl.Region);
    nRegions = numel(regions);

    %% --- Prepare stats per region ---
    means = zeros(nRegions, 1);
    stds  = zeros(nRegions, 1);

    for i = 1:nRegions
        values = tbl.InhibitionOnsetTime(strcmp(tbl.Region, regions{i}));
        means(i) = mean(values);
        stds(i) = std(values);
    end

    %% --- Create Figure ---
    figure('Name', 'Inhibition Summary Plots', 'Color', 'w', 'Position', [100 100 1000 400]);

    % --- Bar Plot with Error Bars ---
    subplot(1,2,1)
    bar(means, 'FaceColor', [0.2 0.4 0.6]);
    hold on;
    errorbar(1:nRegions, means, stds, 'k.', 'LineWidth', 1.5);
    hold off;
    set(gca, 'XTick', 1:nRegions, 'XTickLabel', regions, 'XTickLabelRotation', 45, 'FontSize', 11);
    ylabel('Mean Inhibition Onset Time (s)', 'FontSize', 12);
    title('Bar Plot with Error Bars', 'FontSize', 14);
    grid on;

    % --- Boxplot ---
    subplot(1,2,2)
    boxplot(tbl.InhibitionOnsetTime, tbl.Region);
    ylabel('Inhibition Onset Time (s)');
    title('Boxplot by Brain Region');
    xtickangle(45);
    set(gca, 'FontSize', 11);

    %% --- Save figure ---
    saveas(gcf, 'inhibition_summary_bar_box.png');
    disp('📊 Saved figure as inhibition_summary_bar_box.png');
end
