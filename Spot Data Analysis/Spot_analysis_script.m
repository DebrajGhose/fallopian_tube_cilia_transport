%% plot_msd_from_tracks_refined_fit.m
% Modified to plot all conditions on the same graph for each bead size
% Updated: Uses Non-Linear Least Squares for a visually tighter fit
% V2: Added D values, theoretical MSD, and axis square
clear all
close all

% --- USER-DEFINED PARAMETERS ---
ylim_5um   = [0 35]; 
ylim_100um = [0 2.1]; 
xlim_5um   = [0 300];
xlim_100um = [0 300];
deltaT = 3.0;  % Assumed time between frames in seconds

% NEW PARAMETERS FOR FILTERING
skip_first_lag = true;  % Skip lag=1 to avoid initial tracking artifacts
max_displacement_per_frame = 20;  % Maximum allowed displacement in microns per frame

% PARAMETERS FOR DIFFUSION ANALYSIS
fit_fraction = 0.8;  % SUGGESTED: Keep this low (0.25) to fit the initial diffusion before drift/confinement
fit_max_time = Inf;   % Time limit
show_power_law_fit = true; 
extend_fit_line = true;

% --- THEORETICAL DIFFUSION PARAMETERS ---
show_theoretical = true;  % Set to false to hide theoretical lines
T_kelvin = 298;           % Temperature in Kelvin (25°C)
eta = 8.9e-4;             % Dynamic viscosity of water at 25°C (Pa·s)
kB = 1.380649e-23;        % Boltzmann constant (J/K)
% --------------------------------

rootDir = pwd;
beadFolders = {'5 um redo tracks','100 um tracks'};
conditions  = {'Empty','No Hormones','2 nM E2'};

% Define colors and styles
conditionColors = {'b', 'r', 'g'};
conditionStyles = {'-', '-', '-'};
conditionMarkers = {'none', 'none', 'none'};

figure('Color','w', 'Position', [100, 100, 1200, 500]); 
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Track statistics
total_rejected = 0;
total_accepted = 0;

% Store diffusion analysis results
diffusion_results = struct();
result_counter = 1;

for b = 1:numel(beadFolders)
    nexttile(b); hold on; grid off;
    
    beadPath = findSubdir(rootDir, beadFolders{b});
    if isempty(beadPath)
        warning('Missing bead-size folder: %s', beadFolders{b});
        continue;
    end
    
    plotHandles = [];
    legendLabels = {};
    
    for c = 1:numel(conditions)
        condPath = findSubdir(beadPath, conditions{c});
        if isempty(condPath), continue; end
        
        csvs = dir(fullfile(condPath, '*.csv'));
        if isempty(csvs), continue; end
        
        allD2 = {};
        rejected_count = 0;
        accepted_count = 0;
        
        for k = 1:numel(csvs)
            filePath = fullfile(csvs(k).folder, csvs(k).name);
            T = readTrackMateCSV(filePath);
            
            if ~ismember('FRAME', T.Properties.VariableNames)
                T.FRAME = round(T.POSITION_T / deltaT);
            end

            T = T(isfinite(T.TRACK_ID) & isfinite(T.POSITION_X) & ...
                  isfinite(T.POSITION_Y) & isfinite(T.FRAME), :);
            if height(T) < 2, continue; end
            
            trackIDs = unique(T.TRACK_ID);
            for ti = 1:length(trackIDs)
                idx = (T.TRACK_ID == trackIDs(ti));
                Ttrack = T(idx, :);
                if height(Ttrack) < 2, continue; end
                Ttrack = sortrows(Ttrack, 'FRAME');
                
                x = Ttrack.POSITION_X;
                y = Ttrack.POSITION_Y;
                frames = Ttrack.FRAME;
                n = height(Ttrack);
                
                for i = 1:(n-1)
                    for j = (i+1):n
                        lag = frames(j) - frames(i);
                        if skip_first_lag && lag == 1, continue; end
                        
                        dx = x(j) - x(i);
                        dy = y(j) - y(i);
                        displacement = sqrt(dx^2 + dy^2);
                        
                        if (displacement / lag) > max_displacement_per_frame
                            rejected_count = rejected_count + 1;
                            continue;
                        end
                        
                        accepted_count = accepted_count + 1;
                        d2 = dx^2 + dy^2;
                        
                        if lag > 0
                            if lag > length(allD2)
                                allD2{lag} = d2;
                            else
                                allD2{lag} = [allD2{lag}; d2];
                            end
                        end
                    end
                end
            end
        end
        
        total_rejected = total_rejected + rejected_count;
        total_accepted = total_accepted + accepted_count;
        
        if isempty(allD2), continue; end
        
        maxLag = length(allD2);
        lags = (1:maxLag)';
        msd = zeros(maxLag, 1);
        
        for lag = 1:maxLag
            if ~isempty(allD2{lag})
                msd(lag) = mean(allD2{lag});
            end
        end
        
        valid_idx = msd > 0;
        msd = msd(valid_idx);
        lags = lags(valid_idx);
        time = lags * deltaT;
        
        time = [0; time];
        msd = [0; msd];
        
        h = plot(time, msd, 'LineStyle', conditionStyles{c}, ...
                 'Color', conditionColors{c}, ...
                 'LineWidth', 2, ...
                 'Marker', conditionMarkers{c});
        
        % --- UPDATED FIT ALGORITHM (NON-LINEAR LEAST SQUARES) ---
        D_coeff = NaN;  % Initialize D_coeff
        alpha = NaN;    % Initialize alpha
        
        if length(time) > 3
            non_zero_idx = time > 0;
            time_fit = time(non_zero_idx);
            msd_fit = msd(non_zero_idx);
            
            if fit_max_time < Inf
                mask = time_fit <= fit_max_time;
                time_fit = time_fit(mask);
                msd_fit = msd_fit(mask);
            end
            
            n_fit = max(3, round(length(time_fit) * fit_fraction));
            n_fit = min(n_fit, length(time_fit));
            time_fit = time_fit(1:n_fit);
            msd_fit = msd_fit(1:n_fit);
            
            if length(time_fit) >= 3
                % Step 1: Initial Guess using Log-Log (Standard Polyfit)
                log_t = log10(time_fit);
                log_msd = log10(msd_fit);
                p_guess = polyfit(log_t, log_msd, 1);
                alpha_guess = p_guess(1);
                D_guess = 10^p_guess(2);
                
                % Step 2: Refine using fminsearch (Non-linear minimization)
                % This minimizes sum of squared errors in LINEAR space
                % Model: MSD = D * t^alpha
                modelFun = @(b,t) b(1) .* t.^b(2); 
                % Cost Function: Sum of squared residuals
                costFun = @(b) sum((msd_fit - modelFun(b, time_fit)).^2);
                
                options = optimset('Display','off', 'TolX', 1e-6);
                b_opt = fminsearch(costFun, [D_guess, alpha_guess], options);
                
                D_coeff = b_opt(1);
                alpha = b_opt(2);
                
                % Calculate R-squared
                msd_predicted = D_coeff * time_fit.^alpha;
                ss_res = sum((msd_fit - msd_predicted).^2);
                ss_tot = sum((msd_fit - mean(msd_fit)).^2);
                r_squared = 1 - (ss_res / ss_tot);
                
                diffusion_results(result_counter).bead_size = beadFolders{b};
                diffusion_results(result_counter).condition = conditions{c};
                diffusion_results(result_counter).alpha = alpha;
                diffusion_results(result_counter).D_coeff = D_coeff;
                diffusion_results(result_counter).r_squared = r_squared;
                diffusion_results(result_counter).time_range_fit = [time_fit(1), time_fit(end)];
                result_counter = result_counter + 1;
                
                if show_power_law_fit
                    if extend_fit_line
                        fit_time_start = time_fit(1);
                        fit_time_end = max(time_fit(end), time(end));
                        fit_time = logspace(log10(fit_time_start), log10(fit_time_end), 100);
                    else
                        fit_time = logspace(log10(time_fit(1)), log10(time_fit(end)), 50);
                    end
                    fit_msd = D_coeff * fit_time.^alpha;
                    
                    % Dotted black line
                    plot(fit_time, fit_msd, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                end
            end
        end
        
        % --- UPDATED LEGEND WITH D AND ALPHA ---
        if ~isnan(alpha) && ~isnan(D_coeff)
            legendLabels{end+1} = sprintf('%s (D=%.2e, α=%.2f)', conditions{c}, D_coeff, alpha);
        else
            legendLabels{end+1} = sprintf('%s (n/a)', conditions{c});
        end
        
        plotHandles = [plotHandles, h];
    end
    
    % --- ADD THEORETICAL MSD LINE ---
    if show_theoretical
        % Determine bead diameter based on folder name
        if contains(beadFolders{b}, '5 um', 'IgnoreCase', true)
            diameter_um = 5;    % 5 micron diameter
            time_theory = linspace(0, xlim_5um(2), 200);
        elseif contains(beadFolders{b}, '100 um', 'IgnoreCase', true)
            diameter_um = 100;  % 100 micron diameter
            time_theory = linspace(0, xlim_100um(2), 200);
        else
            diameter_um = NaN;
            time_theory = [];
        end
        
        if ~isnan(diameter_um) && ~isempty(time_theory)
            % Calculate theoretical diffusion coefficient using Stokes-Einstein
            radius_m = (diameter_um / 2) * 1e-6;  % Convert radius to meters
            D_theory_m2s = kB * T_kelvin / (6 * pi * eta * radius_m);  % D in m²/s
            D_theory_um2s = D_theory_m2s * 1e12;  % Convert to μm²/s
            
            % Calculate theoretical MSD for 2D diffusion: MSD = 4*D*t
            msd_theory = 4 * D_theory_um2s * time_theory;
            
            % Plot theoretical line
            h_theory = plot(time_theory, msd_theory, 'm--', 'LineWidth', 2);
            plotHandles = [plotHandles, h_theory];
            legendLabels{end+1} = sprintf('Theory (D=%.2e μm²/s)', D_theory_um2s);
        end
    end
    
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('MSD (\mum^2)', 'Interpreter','tex', 'FontSize', 12);
    
    if contains(beadFolders{b}, '5 um', 'IgnoreCase', true)
        title('5 μm Beads', 'FontSize', 14, 'FontWeight', 'bold');
        xlim(xlim_5um); ylim(ylim_5um);
    elseif contains(beadFolders{b}, '100 um', 'IgnoreCase', true)
        title('100 μm Beads', 'FontSize', 14, 'FontWeight', 'bold');
        xlim(xlim_100um); ylim(ylim_100um);
    else
        title(beadFolders{b}, 'Interpreter','none', 'FontSize', 14);
    end
    
    % --- AXIS SQUARE ---
    axis square;
    
    if ~isempty(plotHandles)
        legend(plotHandles, legendLabels, 'Location', 'northwest', 'FontSize', 9, 'Box', 'on', 'Interpreter', 'none');
    end
end

% Print results
fprintf('\n=== DIFFUSION EXPONENT ANALYSIS (Linear Least Squares) ===\n');
fprintf('%-20s %-15s %-8s %-12s %-8s\n', 'Bead Size', 'Condition', 'α', 'D', 'R²');
fprintf('------------------------------------------------------------------\n');
for i = 1:length(diffusion_results)
    if isfield(diffusion_results, 'alpha') && i <= numel(fieldnames(diffusion_results))
        fprintf('%-20s %-15s %-8.3f %-12.3e %-8.3f\n', ...
            diffusion_results(i).bead_size, diffusion_results(i).condition, ...
            diffusion_results(i).alpha, diffusion_results(i).D_coeff, diffusion_results(i).r_squared);
    end
end

% Print theoretical values
fprintf('\n=== THEORETICAL DIFFUSION COEFFICIENTS (Stokes-Einstein, T=%.0f K) ===\n', T_kelvin);
fprintf('5 μm sphere in water:   D = %.3e μm²/s\n', kB * T_kelvin / (6 * pi * eta * 2.5e-6) * 1e12);
fprintf('100 μm sphere in water: D = %.3e μm²/s\n', kB * T_kelvin / (6 * pi * eta * 50e-6) * 1e12);
fprintf('\n');

sgtitle('Mean Squared Displacement Analysis', 'FontSize', 16, 'FontWeight', 'bold');

%% Helper functions
function p = findSubdir(parent, targetName)
    p = ""; if ~isfolder(parent), return; end
    L = dir(parent); isdir = [L.isdir]; names = {L(isdir).name};
    match = strcmpi(names, targetName);
    if any(match), p = fullfile(parent, names{find(match,1)}); end
end

function T = readTrackMateCSV(filePath)
    opts = detectImportOptions(filePath, 'NumHeaderLines', 0);
    opts.VariableNamesLine = 1;
    opts.DataLines = [4, inf];
    enforce = {'TRACK_ID','POSITION_X','POSITION_Y','POSITION_T'};
    for i = 1:numel(enforce)
        if any(strcmp(opts.VariableNames, enforce{i})), opts = setvartype(opts, enforce{i}, 'double'); end
    end
    T = readtable(filePath, opts);
end