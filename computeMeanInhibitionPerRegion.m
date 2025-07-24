function regionMeans = computeMeanInhibitionPerRegion()
% computeMeanInhibitionPerRegion()
% Loads 'inhibition_onset_summary.csv' and computes the mean inhibition time per region.
% Returns a table with:
% - Region
% - MeanInhibitionTime (in seconds)
% Also saves it as 'mean_inhibition_per_region.csv'

    % Load the CSV
    fname = 'inhibition_onset_summary.csv';
    if ~isfile(fname)
        error('File "%s" not found.', fname);
    end

    tbl = readtable(fname);

    % Remove NaNs if any got in (just in case)
    validRows = ~isnan(tbl.InhibitionOnsetTime);
    tbl = tbl(validRows, :);

    % Group by region and compute mean
    regionNames = unique(tbl.Region);
    meanTimes = [];

    for i = 1:numel(regionNames)
        region = regionNames{i};
        times = tbl.InhibitionOnsetTime(strcmp(tbl.Region, region));
        meanTime = mean(times);
        meanTimes = [meanTimes; {region, meanTime}];
    end

    % Create table
    regionMeans = cell2table(meanTimes, ...
        'VariableNames', {'Region', 'MeanInhibitionTime'});

    % Save to CSV
    writetable(regionMeans, 'mean_inhibition_per_region.csv');
    disp('✅ Saved region means to mean_inhibition_per_region.csv');
end