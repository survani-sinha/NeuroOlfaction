function browsePSTHsForAllVials(dataInput)
% browsePSTHsForAllVials(dataInput)
%   Shows PSTHs (mean ± SEM) for all vials, 16 per page (4×4 layout).
%   Browse with ← and → arrow keys.
%
%   Example usage:
%       browsePSTHsForAllVials('NEW_aPCx_*')
%       browsePSTHsForAllVials(NEW_aPCx_15)

    PST = [-4 8]; binWidth = 0.1;
    edges = PST(1):binWidth:PST(2);
    centers = edges(1:end-1) + binWidth/2;

    allVialRates = containers.Map('KeyType','int32','ValueType','any');

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

    % Aggregate trial-wise PSTHs for each vial
    for d = 1:length(dataList)
        data = dataList{d};
        for u = 1:numel(data)
            for vialIdx = 1:numel(data(u).vials)
                trials = data(u).vials(vialIdx).trials;
                for t = 1:numel(trials)
                    st = trials{t}(:);
                    st = st(st >= PST(1) & st <= PST(2));
                    rates = histcounts(st, edges) / binWidth;
                    if isKey(allVialRates, vialIdx)
                        tmp = allVialRates(vialIdx);
                        tmp{end+1} = rates;
                        allVialRates(vialIdx) = tmp;
                    else
                        allVialRates(vialIdx) = {rates};
                    end
                end
            end
        end
    end

    % Sort vial keys and chunk into groups of 16
    vialKeys = sort(cell2mat(keys(allVialRates)));
    pages = ceil(length(vialKeys)/16);
    currentPage = 1;

    hFig = figure('Name','PSTH Browser','KeyPressFcn',@onKey);
    drawPage();

    function onKey(~, evt)
        switch evt.Key
            case 'rightarrow'
                currentPage = min(currentPage + 1, pages);
            case 'leftarrow'
                currentPage = max(currentPage - 1, 1);
            otherwise
                return;
        end
        drawPage();
    end

    function drawPage()
        clf(hFig);
        tlo = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
        firstIdx = (currentPage - 1)*16 + 1;
        lastIdx = min(currentPage*16, length(vialKeys));
        sgtitle(tlo, sprintf('Average PSTHs – Page %d/%d', currentPage, pages));
        
        for i = firstIdx:lastIdx
            vial = vialKeys(i);
            trialRatesCell = allVialRates(vial);
            trialRates = vertcat(trialRatesCell{:});
    
            meanRate = mean(trialRates, 1);
            semRate = std(trialRates, 0, 1) / sqrt(size(trialRates, 1));
    
            ax = nexttile(tlo);
            axes(ax);  % set current axes
            hold on;
            lineProps = {'Color', 'k'};
            shadedErrorBar(centers, meanRate, semRate, 'lineProps', lineProps);
            xline(0, 'r--');
            xline(1, 'r--');
            title(sprintf('Vial %d', vial));
            xlim(PST);
            ylim([0 max(meanRate + semRate)*1.2]);
        end
    end


end
