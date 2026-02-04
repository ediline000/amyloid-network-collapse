function plot_raster_from_csv(data_folder, sort_by)
time = csvread(fullfile(data_folder, 'time_ms.csv'));
pop_rate = readtable(fullfile(data_folder, 'population_rate.csv'));

% Charger les événements de spike (seulement les vrais spikes!)
events = readtable(fullfile(data_folder, sprintf('spike_events_%s.csv', sort_by)));
if height(events) > 1
    time_spike = events.time_ms;
    rank_spike = events.neuron_rank;
else
    time_spike = []; rank_spike = [];
end

figure('Position', [100, 100, 1000, 400]);

% Raster
subplot(3,1,1);
if ~isempty(time_spike)
    scatter(time_spike, rank_spike, 1, 'k', 'filled');
end
ylabel(['Neuron (sorted by ' upper(sort_by) ')']);
xlim([time(1), time(end)]);
ylim([0,10000]);%min(rank_spike, [], 'omitnan')-1, max(rank_spike, [], 'omitnan')+1]);
set(gca, 'XTickLabel', []);
grid on;

% Population rate lissé
subplot(3,1,2);
plot(pop_rate.time_ms, pop_rate.pop_rate_smooth, 'Color', [0 0.5 0], 'LineWidth', 1);
fill([pop_rate.time_ms; flipud(pop_rate.time_ms)], ...
     [pop_rate.pop_rate_smooth; zeros(size(pop_rate.pop_rate_smooth))], ...
     [0.8 1 0.8], 'EdgeColor', 'none');
ylabel('Smoothed Pop. Rate');
set(gca, 'XTickLabel', []);
%grid on;

% Population rate brut
subplot(3,1,3);
plot(pop_rate.time_ms, pop_rate.pop_rate_raw, 'Color', [0.7 0 0], 'LineWidth', 0.8);
fill([pop_rate.time_ms; flipud(pop_rate.time_ms)], ...
     [pop_rate.pop_rate_raw; zeros(size(pop_rate.pop_rate_raw))], ...
     [1 0.8 0.8], 'EdgeColor', 'none');
ylabel('Raw Pop. Rate');
xlabel('Time (ms)');
grid on;
end
