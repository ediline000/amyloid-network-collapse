function plot_summary_metrics()
% Plot summary metrics (e.g., firing rate, Fano factor) vs disease time
% Assumes folders: spiking_results_t0010, spiking_results_t0100, etc.
% and files: spiking_results_tXXXX/spike_metrics.csv

% --- 1. Define snapshot times ---
snapshot_times = [10, 100, 200, 260, 275, 300, 350];
n_times = numel(snapshot_times);

% --- 2. Define metrics to load and plot ---
metric_names = {
    'mean_firing_rate_all', ...
    'mean_firing_rate_connected', ...
    'connected_spiking', ...
    'population_synchrony', ...
    'fano_factor'
};
metric_titles = {
    'Mean Firing Rate (Hz)', ...
    'Connected Firing Rate (Hz)', ...
    'Spiking Connected Neurons', ...
    'Population Synchrony', ...
    'Fano Factor'
};
colors = [
    0 0 0;      % black
    0 0 1;      % blue
    0 0.5 0;    % green (for connected spiking)
    0 0.8 0;    % lighter green
    1 0.5 0     % orange
];

% --- 3. Load data ---
n_metrics = numel(metric_names);
metric_data = NaN(n_times, n_metrics);

for i = 1:n_times
    t = snapshot_times(i);
    folder = sprintf('spiking_results_t%04d', t);
    file = fullfile(folder, 'spike_metrics.csv');
    
    if exist(file, 'file')
        try
            data = readtable(file);
            % Assume first row contains the scalar metrics
            for j = 1:n_metrics
                if ismember(metric_names{j}, data.Properties.VariableNames)
                    metric_data(i, j) = data.(metric_names{j})(1);
                end
            end
            fprintf('Loaded metrics for t = %d\n', t);
        catch
            warning('Failed to read %s', file);
        end
    else
        fprintf('⚠️ Missing metrics for t = %d\n', t);
    end
end

% --- 4. Plot ---
n_cols = 2;
n_rows = ceil(n_metrics / n_cols);
figure('Position', [100, 100, 1200, 400 * n_rows]);

for j = 1:n_metrics
    subplot(n_rows, n_cols, j);
    y = metric_data(:, j);
    valid = ~isnan(y);
    
    if any(valid)
        plot(snapshot_times(valid), y(valid), 'o-', 'Color', colors(j,:), ...
             'LineWidth', 2, 'MarkerSize', 6);
        title(metric_titles{j}, 'FontWeight', 'bold');
        xlabel('Disease Time (months)');
        ylabel(metric_titles{j});
        grid on;
        box on;
        
        % Optional: vertical line at t=200
        xline(200, '--', 'Color', [0.5 0.5 0.5], 'Alpha', 0.5);
    else
        title(metric_titles{j}, 'FontWeight', 'bold');
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 12);
    end
end

sgtitle('Evolution of Network Spiking Dynamics During Disease Progression', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Save figure
% saveas(gcf, 'spiking_metrics_evolution_MATLAB.png');
% disp('✅ Figure saved as spiking_metrics_evolution_MATLAB.png');
end