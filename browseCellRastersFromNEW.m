function browseCellRastersFromNEW(data)
% browseCellRastersFromNEW(data)
%   Pages through rasters + PSTHs (←/→) for all cells in the data.
%   Layout: 4×4 grid, one panel per vial.

    PST       = [-4 8];
    binWidth  = 0.1;
    odorDur   = 1.0;
    nPanels   = 16;
    maxTrials = 10;

    nUnits = numel(data);
    unitIdx = 1;

    hFig = figure('Name','Raster+PSTH Browser','NumberTitle','off', ...
                  'KeyPressFcn',@onKey);
    redraw();

    function onKey(~,evt)
        switch evt.Key
            case 'leftarrow'
                unitIdx = unitIdx - 1;
            case 'rightarrow'
                unitIdx = unitIdx + 1;
            otherwise
                return;
        end
        unitIdx = mod(unitIdx-1, nUnits) + 1;
        redraw();
    end

    function redraw()
        clf(hFig);
        tlo = tiledlayout(4,4, ...
            'TileSpacing','tight', ...
            'Padding','compact');  % ← tighter fit, less overlap
        sgtitle(tlo, sprintf('Unit %d/%d  (←/→)', unitIdx, nUnits));

        vials = data(unitIdx).vials;
        for vial = 1:nPanels
            axT = nexttile(tlo, vial);
            pos = axT.Position;
            delete(axT);

            % Adjust inner padding to shrink plots slightly
            rasterHeight = pos(4)*0.45;
            psthHeight   = pos(4)*0.4;
            gap          = pos(4)*0.05;

            if vial <= numel(vials)
                trials = vials(vial).trials;
            else
                trials = {};
            end
            useTrials = min(numel(trials), maxTrials);
            if useTrials == 0, continue; end

            % Raster plot (top)
            axR = axes('Position', [pos(1), pos(2)+psthHeight+gap, pos(3), rasterHeight]);
            hold(axR, 'on');
            for t = 1:useTrials
                st = trials{t}(:);
                st = st(st >= PST(1) & st <= PST(2));
                plot(axR, st, t*ones(size(st)), 'k.', 'MarkerSize', 3);
            end
            xline(axR, [0 odorDur], 'r-');
            xlim(axR, PST);
            ylim(axR, [0.5 useTrials+0.5]);
            axis(axR, 'off');
            title(axR, sprintf('Vial %d', vial), 'FontSize', 8);

            % PSTH plot (bottom)
            axP = axes('Position', [pos(1), pos(2), pos(3), psthHeight]);
            allSpikes = [];
            for t = 1:useTrials
                allSpikes = [allSpikes; trials{t}(:)]; %#ok<AGROW>
            end
            if ~isempty(allSpikes)
                edges   = PST(1):binWidth:PST(2);
                counts  = histcounts(allSpikes, edges);
                rates   = counts / (useTrials * binWidth);
                centers = edges(1:end-1) + binWidth/2;
                bar(axP, centers, rates, 'hist');
            end
            xline(axP, [0 odorDur], 'r-');
            xlim(axP, PST);
            axis(axP, 'on');
            set(axP, 'FontSize', 6);
            xlabel(axP, 'Time (s)', 'FontSize', 6);
            ylabel(axP, 'Hz', 'FontSize', 6);
        end
    end
end
