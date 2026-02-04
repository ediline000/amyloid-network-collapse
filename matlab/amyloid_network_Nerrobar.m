function amyloid_network_Nerrobar
    clc
    %% Network Parameters
    N_values =[100,8000]; % Different number of neurons to test
    k = 20;                % Average degree
    p_rewire = 0.1;        % Rewiring probability
    t_max = 500;           % Simulation time (months)
    dt = 1;                % Time step
    n_realizations = 5;    % Number of trials
    
    %% Biological Parameters (same as before)
    A_max = 150;                % Maximum amyloid level
    k_A = 0.02;                 % Amyloid growth rate
    A_th = 25;                  % Amyloid toxicity threshold
    clearance_rate = 0.005;     % Microglial clearance rate
    W_0 = 1;                    % Initial synaptic weight
    W_min = 0.5 * W_0 * exp(-1); % Minimum weight threshold
    tau_synapse = 120;          % Synaptic decay time constant
    base_repair_rate = 0.02;    % Baseline repair probability
    repair_strength = 0.03;     % Weight recovery factor
    
    %% Initialize Storage
    n_steps = t_max/dt;
    time_points = linspace(0,t_max,n_steps);
    
    % Storage for metrics for each N
    amyloid = zeros(n_steps, n_realizations, length(N_values));
    lscc = zeros(n_steps, n_realizations, length(N_values));
    efficiency = zeros(n_steps, n_realizations, length(N_values));
    repair_events = zeros(n_steps, n_realizations, length(N_values));
    V = zeros(n_steps, length(N_values)); % Store V for each N
    
    %% Main Simulation (same as before)
    for n_idx = 1:length(N_values)
        N = N_values(n_idx);
        for realization = 1:n_realizations
            fprintf('N = %d, Realization %d/%d...\n', N, realization, n_realizations);
            
            % Create small-world network
            G = wattsStrogatz(N, k, p_rewire);
            G = make_directed(G, W_0);
            
            % Initialize variables
            A = 4; % Initial amyloid level
            
            for i = 1:n_steps
                %% Amyloid Dynamics (with microglial clearance)
                A = A + dt * (k_A * A * (1 - A / A_max));
                A = A - dt * clearance_rate * A;
                A = max(0, A);
                amyloid(i, realization, n_idx) = A;
                
                %% Synaptic Degradation
                p_toxicity = 1 - exp(-0.001 * max(0, A - A_th));
                edges_to_remove = [];
                for e = 1:numedges(G)
                    if rand() < p_toxicity
                        G.Edges.Weight(e) = W_0 * exp(-(i * dt) / tau_synapse);
                    end
                    if G.Edges.Weight(e) < W_min
                        edges_to_remove = [edges_to_remove; e];
                    end
                end
                G = rmedge(G, edges_to_remove);
                
                %% Network Repair Mechanisms
                current_repair_rate = base_repair_rate * (1 - exp(-0.01 * (A_th - A))) / sqrt(N);
                [G, n_repairs] = repair_network(G, W_0, current_repair_rate, repair_strength);
                repair_events(i, realization, n_idx) = n_repairs;
                
                %% Compute Network Metrics
                [~, sizes] = conncomp(G);
                lscc(i, realization, n_idx) = max(sizes) / N;
                
                % Global efficiency
                adj_mat = full(adjacency(G, 'weighted'));
                efficiency(i, realization, n_idx) = global_efficiency_fun(adj_mat);
            end
        end
        
        %% Calculate V for each N
        mean_lscc = mean(lscc(:,:,n_idx), 2);
        mean_sc = mean(lscc(:,:,n_idx).^2, 2);
        V(:,n_idx) = (mean_sc ./ mean_lscc) - mean_lscc;
        
        % Find critical time
        [~, tc_idx] = max(V(:,n_idx));
        TC = time_points(tc_idx);
        fprintf('Critical time T_c for N = %d: %.2f\n', N, TC);
    end
    
    %% Calculate Means and Standard Deviations
    mean_amyloid = mean(amyloid, 2);
    std_amyloid = std(amyloid, 0, 2);
    mean_efficiency = mean(efficiency, 2);
    std_efficiency = std(efficiency, 0, 2);
    mean_repairs = mean(repair_events, 2);
    std_repairs = std(repair_events, 0, 2);
    mean_lscc = mean(lscc, 2);
    std_lscc = std(lscc, 0, 2);
    
    %% Plot all N values on same graphs with error bars
    colors = lines(length(N_values)); % Different colors for each N
    
    % Figure 1: Amyloid Levels
    figure(7);
    hold on;
    for n_idx = 1:length(N_values)
        errorbar(time_points, squeeze(mean_amyloid(:,:,n_idx)), squeeze(std_amyloid(:,:,n_idx)), ...
                'Color', colors(n_idx,:), 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(N_values(n_idx))]);
    end
    title('Amyloid Dynamics');
    xlabel('Time (months)');
    ylabel('Amyloid Level');
    legend('show');
    grid on;
    hold off;
    
    % Figure 2: LSCC Fraction
    figure(8);
    hold on;
    for n_idx = 1:length(N_values)
        errorbar(time_points, squeeze(mean_lscc(:,:,n_idx)), squeeze(std_lscc(:,:,n_idx)), ...
                'Color', colors(n_idx,:), 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(N_values(n_idx))]);
    end
    title('Network Connectivity (LSCC Fraction)');
    xlabel('Time (months)');
    ylabel('LSCC Fraction');
    legend('show');
    grid on;
    hold off;
    
    % Figure 3: Global Efficiency
    figure(9);
    hold on;
    for n_idx = 1:length(N_values)
        errorbar(time_points, squeeze(mean_efficiency(:,:,n_idx)), squeeze(std_efficiency(:,:,n_idx)), ...
                'Color', colors(n_idx,:), 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(N_values(n_idx))]);
    end
    title('Network Efficiency');
    xlabel('Time (months)');
    ylabel('Global Efficiency');
    legend('show');
    grid on;
    hold off;
    
    % Figure 4: Repair Events
    figure(10);
    hold on;
    for n_idx = 1:length(N_values)
        errorbar(time_points, squeeze(mean_repairs(:,:,n_idx)), squeeze(std_repairs(:,:,n_idx)), ...
                'Color', colors(n_idx,:), 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(N_values(n_idx))]);
    end
    title('Network Repair Activity');
    xlabel('Time (months)');
    ylabel('Repair Events');
    legend('show');
    grid on;
    hold off;
    
    % Figure 5: V values
    figure(11);
    hold on;
    for n_idx = 1:length(N_values)
        plot(time_points, V(:,n_idx), 'Color', colors(n_idx,:), 'LineWidth', 1.5, 'DisplayName', ['N = ' num2str(N_values(n_idx))]);
    end
    title('V Values');
    xlabel('Time (months)');
    ylabel('V');
    legend('show');
    grid on;
    hold off;
end
function G = wattsStrogatz(N, k, beta)
    % Create Watts-Strogatz small-world network
    G = graph();
    G = addnode(G, N);
    
    % Create ring lattice
    for i = 1:N
        for j = 1:k/2
            neighbor = mod(i + j - 1, N) + 1;
            G = addedge(G, i, neighbor, 1);
        end
    end
    
    % Rewire edges with probability beta
    edges = table2array(G.Edges);
    for e = 1:size(edges, 1)
        if rand() < beta
            G = rmedge(G, edges(e,1), edges(e,2));
            new_target = randi(N);
            while new_target == edges(e,1) || findedge(G, edges(e,1), new_target) > 0
                new_target = randi(N);
            end
            G = addedge(G, edges(e,1), new_target, 1);
        end
    end
end

function G = make_directed(G, W_0)
    % Convert undirected graph to directed with random orientation
    [src, tgt] = findedge(G);
    for e = 1:numedges(G)
        if rand() < 0.5
            G = rmedge(G, src(e), tgt(e));
            G = addedge(G, tgt(e), src(e), W_0);
        else
            G.Edges.Weight(e) = W_0;
        end
    end
end

function [G, n_repairs] = repair_network(G, W_0, p_repair, repair_strength)
    % Attempt repair on each existing synapse
    n_repairs = 0;
    for e = 1:numedges(G)
        if rand() < p_repair
            G.Edges.Weight(e) = min(W_0, G.Edges.Weight(e) * (1 + repair_strength));
            n_repairs = n_repairs + 1;
        end
    end
    
    % Add new connections with probability proportional to repair rate
    if rand() < p_repair/2
        nodes = randperm(numnodes(G), 2);
        if ~findedge(G, nodes(1), nodes(2))
            G = addedge(G, nodes(1), nodes(2), W_0*0.5); % New weak synapse
            n_repairs = n_repairs + 1;
        end
    end
end
function lscc_ratio = compute_lscc(G, N)
    if isempty(G)
        lscc_ratio = 0;
        return;
    end
    
    % Create a directed graph object
    G_directed = digraph(G);

    % Find strongly connected components
    [bins, ~] = conncomp(G_directed, 'Type', 'strong');

    % Count the sizes of each component
    component_sizes = histcounts(bins, 'BinMethod', 'integers');
    
    % Find the size of the largest strongly connected component
    largest_scc_size = max(component_sizes);
    
    % Compute the ratio of the largest SCC size to N
    lscc_ratio = largest_scc_size / N;
end
function GE = global_efficiency_fun(G)
     % Replace zeros with infinity for computing efficiency
    G(G == 0) = Inf;
    
    % Create the graph with weights as the inverse of G
    graph = G;
    graph(graph > 0) = 1 ./ graph(graph > 0);
    
    % Create a directed graph object
    G_directed = digraph(graph);
    
    % Compute shortest paths using the distances function
    shortest_paths = distances(G_directed);
    
    % Set the diagonal to infinity
    shortest_paths(1:size(shortest_paths, 1) + 1:end) = Inf;
    
    % Compute efficiency
    efficiency = zeros(size(shortest_paths));
    efficiency(shortest_paths > 0) = 1 ./ shortest_paths(shortest_paths > 0);
    
    % Calculate Global Efficiency
    N = size(G, 1);
    GE = sum(efficiency(:)) / (N * (N - 1));
end