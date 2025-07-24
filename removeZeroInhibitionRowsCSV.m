function removeZeroInhibitionRowsCSV(inputFile)
    % removeZeroInhibitionRowsCSV removes rows with InhibitionOnsetTime == 0
    % from a CSV file and saves the result as cleanInhibition.csv.
    %
    % Input:
    %   inputFile - string, name of the input CSV file (e.g., 'data.csv')

    % Read the CSV file into a table
    data = readtable(inputFile);

    % Check that the required column exists
    if ~ismember('InhibitionOnsetTime', data.Properties.VariableNames)
        error('Column "InhibitionOnsetTime" not found in the input file.');
    end

    % Filter the rows
    filteredData = data(data.InhibitionOnsetTime ~= 0, :);

    % Save the filtered table to a new CSV file
    writetable(filteredData, 'cleanInhibition.csv');

    % Optional: Display confirmation
    disp('Filtered data saved as cleanInhibition.csv');
end
