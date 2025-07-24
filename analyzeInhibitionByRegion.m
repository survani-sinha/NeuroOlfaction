function avgTable = analyzeInhibitionByRegion(csvFile)
    % analyzeInhibitionByRegion computes mean, std, and SEM of InhibitionOnsetTime per Region,
    % then generates bar plot with error bars and boxplot.
    
    % Read table
    data = readtable(csvFile);

    % Validate necessary columns
    if ~all(ismember({'Region', 'InhibitionOnsetTime'}, data.Properties.VariableNames))
        error('CSV must contain columns "Region" and "InhibitionOnsetTime".');
    end

    % Use groupsummary to compute mean and std
    stats = groupsummary(data, 'Region', {'mean', 'std'}, 'InhibitionOnsetTime');

    % Extract values
    regions = stats.Region;
    means   = stats.mean_InhibitionOnsetTime;
    stds    = stats.std_InhibitionOnsetTime;
    ns      = stats.GroupCount;
    sems    = stds ./ sqrt(ns);

    % Create output table
    avgTable = table(regions, means, stds, ns, sems, ...
        'VariableNames', {'Region', 'Mean', 'StdDev', 'N', 'SEM'});

    % Save to CSV
    writetable(avgTable, 'avgInhibitionByRegion_withError.csv');
    disp(avgTable);

    %% --- Plot 1: Bar plot with error bars (±SEM) ---
    figure;
    bar(categorical(avgTable.Region), avgTable.Mean, 'FaceColor', [0.4 0.6 0.85]);
    hold on;
    errorbar(categorical(avgTable.Region), avgTable.Mean, avgTable.SEM, ...
        'k.', 'LineWidth', 1.5, 'CapSize', 10);
    title('Average Inhibition Onset Time by Region (±SEM)');
    xlabel('Region');
    ylabel('Inhibition Onset Time (s)');
    xtickangle(45);
    grid on;
    hold off;

    %% --- Plot 2: Boxplot ---
    figure;
    boxplot(data.InhibitionOnsetTime, data.Region);
    title('Inhibition Onset Time Distribution by Region');
    xlabel('Region');
    ylabel('Inhibition Onset Time (s)');
    grid on;
end
