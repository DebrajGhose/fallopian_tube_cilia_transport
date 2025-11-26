% Mechanism-Based Transport Simulation (Unified Physics + Visualization)
close all; clear; clc;

%% 1. Parameters
dt = 0.05;                  % seconds
simulation_time = 60;       % seconds
n_steps = simulation_time / dt;
n_simulations = 100;         % More particles for smooth averages

% Bead Sizes
bead_configs = [5, 100];    

% Conditions
conditions = struct('name', {'Empty', 'Active'}, ...
                    'has_cilia', {false, true}, ...
                    'color', {'b', 'g'});

% Physics Constants (In MICRON Units)
kB = 1.38e-23;              
T = 298;                    
eta = 1e-3;                 % Pa.s
% Drag coefficient factor (6*pi*eta) converted for micron output
drag_factor = 6 * pi * eta * 1e-6; 

% Mechanism Parameters
% 1. Obstacle Parameters (For 5um)
cilia_spacing = 10;          % Average gap between obstacles (um)
cilia_stem_rad = 0.8;       % Radius of the obstacle (um)

% 2. Active Parameters (For 100um)
% The 'kick' velocity provided by cilia beating
kick_velocity = 2.0;        % um/s (Magnitude of active noise)

%% 2. Generate Landscape
% We generate a fixed field of obstacles
L = 500; % Domain size
num_cilia = round((L/cilia_spacing)^2);
obs_x = (rand(num_cilia,1)-0.5)*L;
obs_y = (rand(num_cilia,1)-0.5)*L;
tree = KDTreeSearcher([obs_x, obs_y]);

%% 3. Main Loop
results = struct();

for b_idx = 1:length(bead_configs)
    diam = bead_configs(b_idx);
    rad = diam / 2;
    
    % Calculate Diffusion Coefficient (Micron^2/s)
    gamma = drag_factor * rad; % Drag in SI
    D_si = (kB * T) / gamma;   % m^2/s
    D = D_si * 1e12;           % um^2/s
    
    fprintf('Bead: %d um, D: %.4f um^2/s\n', diam, D);
    
    % REGIME CHECK:
    is_immersed = diam < cilia_spacing * 2; 

    for c_idx = 1:length(conditions)
        cond = conditions(c_idx);
        msd_sum = zeros(1, n_steps);
        
        for sim = 1:n_simulations
            x = 0; y = 0;
            traj_x = zeros(1, n_steps);
            traj_y = zeros(1, n_steps);
            
            for t = 1:n_steps
                
                % 1. Base Brownian Motion (Always present)
                dx = sqrt(2*D*dt) * randn();
                dy = sqrt(2*D*dt) * randn();
                
                % 2. Cilia Interactions
                if cond.has_cilia
                    if is_immersed
                        % --- REGIME: HINDERED DIFFUSION ---
                        % Proposed new position
                        x_new = x + dx;
                        y_new = y + dy;
                        
                        % Check for collisions with stems
                        search_r = rad + cilia_stem_rad;
                        idx = rangesearch(tree, [x_new, y_new], search_r);
                        
                        if ~isempty(idx{1})
                            % COLLISION DETECTED: Bounce/Reject
                            dx = dx * -0.5; 
                            dy = dy * -0.5;
                        end
                        
                    else
                        % --- REGIME: ACTIVE BATH ---
                        % Bead sits on top, receives random kicks
                        angle = rand() * 2 * pi;
                        dx_kick = cos(angle) * kick_velocity * dt;
                        dy_kick = sin(angle) * kick_velocity * dt;
                        
                        dx = dx + dx_kick;
                        dy = dy + dy_kick;
                    end
                end
                
                x = x + dx;
                y = y + dy;
                
                traj_x(t) = x;
                traj_y(t) = y;
            end
            msd_sum = msd_sum + (traj_x.^2 + traj_y.^2);
        end
        results(b_idx, c_idx).time = (1:n_steps)*dt;
        results(b_idx, c_idx).msd = msd_sum / n_simulations;
    end
end

%% 4. Plotting MSD Results
figure('Position', [100, 400, 1000, 400], 'Color', 'w');
tiledlayout(1,2);

% 5 micron Plot
nexttile; hold on; title('5 \mum Beads');
plot(results(1,1).time, results(1,1).msd, 'b', 'LineWidth', 2); % Empty
plot(results(1,2).time, results(1,2).msd, 'g', 'LineWidth', 2); % Active
xlabel('Time (s)'); ylabel('MSD (\mu m^2)'); grid on; legend('Empty', 'Active');

% 100 micron Plot
nexttile; hold on; title('100 \mum Beads');
plot(results(2,1).time, results(2,1).msd, 'b', 'LineWidth', 2); % Empty
plot(results(2,2).time, results(2,2).msd, 'g', 'LineWidth', 2); % Active
xlabel('Time (s)'); ylabel('MSD (\mu m^2)'); grid on; legend('Empty', 'Active');


%% 5. VISUALIZATION OF SCALE (The requested Figure)
figure('Position', [100, 50, 600, 600], 'Color', 'w');
axis equal; hold on; box on;
title('Scale Model: Cilia Forest vs. Bead Sizes', 'FontSize', 14);
xlabel('Position (\mu m)'); ylabel('Position (\mu m)');

% A. Draw the Cilia Patch (Zoomed in area)
vis_range = 80; % View window size
xlim([-vis_range/2, vis_range/2]);
ylim([-vis_range/2, vis_range/2]);

% Filter cilia to only draw those in view
in_view = abs(obs_x) < vis_range/2 + 5 & abs(obs_y) < vis_range/2 + 5;
cx = obs_x(in_view);
cy = obs_y(in_view);

% Draw Cilia Stems (Gray circles)
for i = 1:length(cx)
    rectangle('Position', [cx(i)-cilia_stem_rad, cy(i)-cilia_stem_rad, ...
              cilia_stem_rad*2, cilia_stem_rad*2], ...
              'Curvature', [1 1], 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
end

% B. Draw 5um Bead (Blue)
% Place it in a gap near the center
r5 = 2.5;
rectangle('Position', [-10-r5, -10-r5, r5*2, r5*2], ...
          'Curvature', [1 1], 'FaceColor', 'b', 'EdgeColor', 'k');
text(-10, -10+r5+2, '5 \mum', 'Color', 'b', 'HorizontalAlignment', 'center');

% C. Draw 100um Bead (Green, transparent)
% Place it offset so it dominates the view
r100 = 50;
rectangle('Position', [20-r100, 20-r100, r100*2, r100*2], ...
          'Curvature', [1 1], 'FaceColor', 'g', 'EdgeColor', 'k', ...
          'FaceAlpha', 0.3); % FaceAlpha requires MATLAB R2014b or later
text(20, 20, '100 \mum', 'Color', [0 0.5 0], 'HorizontalAlignment', 'center', 'FontSize', 12);

% Add Annotation
text(0, -vis_range/2 + 5, 'Gray Dots = Cilia Stems', 'HorizontalAlignment', 'center');