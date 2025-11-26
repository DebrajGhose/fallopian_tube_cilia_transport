%% plot_normalized_traces_from_tracks.m
clear all; close all;

% --- USER PARAMETERS ---
deltaT = 3.0;                 % seconds per frame (only used if FRAME missing)
maxTracesToPlot = 10;         % per panel
beadFolders = {'5 um redo tracks','100 um tracks'};
conditions  = {'Empty','No Hormones','2 nM E2'};
rootDir = pwd;
% ------------------------

figure('Color','w');
tiledlayout(numel(beadFolders), numel(conditions), ...
    'TileSpacing','compact','Padding','compact');

panel = 0;
for b = 1:numel(beadFolders)
    beadPath = findSubdir(rootDir, beadFolders{b});
    for c = 1:numel(conditions)
        panel = panel + 1; nexttile(panel); hold on; grid on;
        title(sprintf('%s — %s', beadFolders{b}, conditions{c}), 'Interpreter','none');
        xlabel('\DeltaX'); ylabel('\DeltaY');

        if isempty(beadPath), axis off; continue; end
        condPath = findSubdir(beadPath, conditions{c});
        if isempty(condPath), axis off; continue; end

        csvs = dir(fullfile(condPath, '*.csv'));
        if isempty(csvs), axis off; continue; end

        % --- Collect all tracks across CSVs ---
        tracks = {};  % each cell: table for one track
        for k = 1:numel(csvs)
            filePath = fullfile(csvs(k).folder, csvs(k).name);
            T = readTrackMateCSV(filePath);

            % FRAME fallback if missing
            if ~ismember('FRAME', T.Properties.VariableNames)
                T.FRAME = round(T.POSITION_T / deltaT);
            end

            T = T(isfinite(T.TRACK_ID) & isfinite(T.POSITION_X) & ...
                  isfinite(T.POSITION_Y) & isfinite(T.FRAME), :);
            if height(T) < 2, continue; end

            trackIDs = unique(T.TRACK_ID);
            for ti = 1:numel(trackIDs)
                idx = (T.TRACK_ID == trackIDs(ti));
                Ttrack = sortrows(T(idx,:), 'FRAME');
                if height(Ttrack) >= 2
                    tracks{end+1} = Ttrack; %#ok<AGROW>
                end
            end
        end

        if isempty(tracks)
            text(0.5,0.5,'No usable tracks','HorizontalAlignment','center');
            axis off; continue;
        end

        % --- Pick up to maxTracesToPlot longest tracks ---
        lens = cellfun(@(Tt) height(Tt), tracks);
        [~, order] = sort(lens, 'descend');
        pick = order(1:min(maxTracesToPlot, numel(order)));

        % --- Plot each selected track, re-anchored at (0,0) ---
        for pidx = pick(:)'
            Tt = tracks{pidx};
            x = Tt.POSITION_X;  y = Tt.POSITION_Y;
            x0 = x(1);          y0 = y(1);
            dx = x - x0;        dy = y - y0;

            plot(dx, dy, '-o', 'MarkerSize', 3);         % full path
            plot(0, 0, 'ks', 'MarkerSize', 5, 'LineWidth', 1);  % origin marker
        end

        axis equal;
        % Optional: auto-tight limits with a small margin
        xl = xlim; yl = ylim;
        pad = 0.05 * max(diff(xl), diff(yl));
        xlim([xl(1)-pad, xl(2)+pad]); ylim([yl(1)-pad, yl(2)+pad]);
        legend(sprintfc('Trace %d', 1:numel(pick)), 'Location','bestoutside'); %#ok<*NOPRT>
    end
end

%% Helper functions (same as before)
function p = findSubdir(parent, targetName)
    p = "";
    if ~isfolder(parent), return; end
    L = dir(parent);
    isdir = [L.isdir];
    names = {L(isdir).name};
    names = names(~ismember(lower(names), {'.','..'}));
    if isempty(names), return; end
    match = strcmpi(names, targetName);
    if any(match)
        p = fullfile(parent, names{find(match,1)});
    else
        match = contains(lower(names), lower(targetName));
        if any(match), p = fullfile(parent, names{find(match,1)}); end
    end
end

function T = readTrackMateCSV(filePath)
    opts = detectImportOptions(filePath, 'NumHeaderLines', 0);
    opts.VariableNamesLine = 1;
    opts.DataLines = [4, inf];  % TrackMate: headers+units in first 3 lines
    enforce = {'TRACK_ID','POSITION_X','POSITION_Y','POSITION_T'};
    for i = 1:numel(enforce)
        if any(strcmp(opts.VariableNames, enforce{i}))
            opts = setvartype(opts, enforce{i}, 'double');
        end
    end
    T = readtable(filePath, opts);
end
