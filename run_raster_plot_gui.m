function plot_population_psth_multiStructs_workspace()
% GUI-based version that lets you select multiple OB and PCx structs from workspace

    % Get all struct variables from workspace
    vars = evalin('base', 'whos');
    structVars = vars(strcmp({vars.class}, 'struct'));

    if isempty(structVars)
        error('No struct variables found in the workspace.');
    end

    % Select OB structs
    ob_idx = listdlg('ListString', {structVars.name}, ...
                     'PromptString', 'Select OB struct(s)', ...
                     'SelectionMode', 'multiple', ...
                     'ListSize', [300 300]);
    if isempty(ob_idx)
        error('No OB struct selected.');
    end
    ob_structs = cellfun(@(n) evalin('base', n), ...
                         {structVars(ob_idx).name}, ...
                         'UniformOutput', false);

    % Select PCx structs
    pcx_idx = listdlg('ListString', {structVars.name}, ...
                      'PromptString', 'Select PCx struct(s)', ...
                      'SelectionMode', 'multiple', ...
                      'ListSize', [300 300]);
    if isempty(pcx_idx)
        error('No PCx struct selected.');
    end
    pcx_structs = cellfun(@(n) evalin('base', n), ...
                          {structVars(pcx_idx).name}, ...
                          'UniformOutput', false);

    % Call main plotting function
    plot_population_psth_from_multiple_structs(ob_structs, pcx_structs);
end
