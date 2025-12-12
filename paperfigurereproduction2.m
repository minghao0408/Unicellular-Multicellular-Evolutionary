%% Fixed Combined Evolutionary Branching and Spatial Distribution Figure
% Corrected layout to show all time labels properly

clear; clc;

%% Load data (assuming already loaded)
if ~exist('adhspec_data', 'var') || ~exist('output_data', 'var')
    fprintf('Loading data files...\n');
    adhspec_data = readmatrix('adhspecfile.txt');
    output_data = readmatrix('outputfile.txt');
end

rows = 200;
cols = 200;
total_frames = size(adhspec_data, 1) / rows;

%% Create figure with corrected layout
figure('Name', 'Fixed: Evolutionary Branching with Spatial Snapshots', 'Position', [100, 100, 1200, 650]);

%% Left panel: Evolutionary branching plot
subplot('Position', [0.08, 0.1, 0.4, 0.8]);  % [left, bottom, width, height]

% Collect data for branching plot
all_times = [];
all_adhesions = [];

% Sample every 3rd frame
sample_frames = 1:3:total_frames;

for frame_idx = 1:length(sample_frames)
    frame = sample_frames(frame_idx);
    start_row = (frame-1) * rows + 1;
    end_row = frame * rows;
    
    current_adhspec = adhspec_data(start_row:end_row, :);
    current_output = output_data(start_row:end_row, :);
    
    cell_positions = current_output == 1;
    if sum(cell_positions(:)) > 10
        adhesion_values = current_adhspec(cell_positions);
        
        % Sample cells
        num_samples = min(30, length(adhesion_values));
        if length(adhesion_values) > num_samples
            sample_idx = randperm(length(adhesion_values), num_samples);
            sampled_adhesions = adhesion_values(sample_idx);
        else
            sampled_adhesions = adhesion_values;
        end
        
        % Calculate actual time
        actual_time = (frame-1) * 1000;
        
        % Store data
        all_times = [all_times; repmat(actual_time, length(sampled_adhesions), 1)];
        all_adhesions = [all_adhesions; sampled_adhesions];
    end
end

% Create scatter plot
scatter(all_adhesions, all_times, 4, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
hold on;

% Add MC/UC boundary
plot([0, 0], [0, 200000], 'r-', 'LineWidth', 2);

% Set proper axes
xlim([-6, 12]);
ylim([0, 200000]);
xlabel('J (S_r)', 'FontSize', 12);
ylabel('Time Steps', 'FontSize', 12);

% Add labels
text(-4, 180000, 'MC', 'FontSize', 14, 'Color', 'blue', 'FontWeight', 'bold');
text(8, 180000, 'UC', 'FontSize', 14, 'Color', 'red', 'FontWeight', 'bold');

% Add φ annotation
text(-5, 10000, 'φ ≈ 90', 'FontSize', 12, 'FontWeight', 'bold');

grid on;
set(gca, 'YDir', 'normal');

%% Right panel: Spatial snapshots with FIXED positions
% Create 3 separate subplot areas with explicit positioning

% Define snapshot parameters
snapshot_times = [50000, 120000, 200000];
snapshot_frames = round(snapshot_times / 1000) + 1;

% Fixed positions for the 3 snapshots [left, bottom, width, height]
positions = [
    0.58, 0.65, 0.25, 0.22;  % Top snapshot
    0.58, 0.38, 0.25, 0.22;  % Middle snapshot  
    0.58, 0.11, 0.25, 0.22   % Bottom snapshot
];

for i = 1:3
    frame = snapshot_frames(i);
    if frame <= total_frames
        start_row = (frame-1) * rows + 1;
        end_row = frame * rows;
        
        current_adhspec = adhspec_data(start_row:end_row, :);
        current_output = output_data(start_row:end_row, :);
        
        % Create spatial image
        spatial_image = zeros(rows, cols);
        for r = 1:rows
            for c = 1:cols
                if current_output(r, c) == 1
                    spatial_image(r, c) = current_adhspec(r, c);
                else
                    spatial_image(r, c) = NaN;
                end
            end
        end
        
        % Create subplot with explicit position
        ax = subplot('Position', positions(i, :));
        
        % Display spatial pattern
        imagesc(spatial_image);
        axis equal; axis tight; axis off;
        
        % Set colormap and limits
        colormap(ax, jet);
        caxis([-20, 10]);
        
        % Add time label - this should now be visible
        title(sprintf('t = %d', snapshot_times(i)), 'FontSize', 11, 'FontWeight', 'bold');
    end
end

%% Add single colorbar for all spatial plots
colorbar_ax = axes('Position', [0.86, 0.2, 0.02, 0.6]);
cb = colorbar(colorbar_ax);
caxis([-20, 10]);
colormap(colorbar_ax, jet);
set(colorbar_ax, 'Visible', 'off');
ylabel(cb, 'Adhesion Value J', 'FontSize', 11);

fprintf('Created fixed layout figure\n');
fprintf('All time labels should now be visible\n');