function amyloid_network_errobars
    clc
    %% Network Parameters
    N = 8000;               % Number of neurons
    k = 20;                % Average degree
    p_rewire = 0.1;        % Rewiring probability
    t_max = 500;           % Simulation time (months)
    dt = 1;                % Time step
    n_realizations = 2;   % Number of trials
    
    %% Biological Parameters
    A_max = 150;           % Maximum amyloid level
    k_A = 0.02;            % Amyloid growth rate
    A_th = 25;             % Amyloid toxicity threshold
    clearance_rate = 0.005; % Microglial clearance rate
    W_0 = 1;               % Initial synaptic weight
    W_min = 0.5 * W_0 * exp(-1); % Minimum weight threshold
    tau_synapse = 120;     % Synaptic decay time constant
    base_repair_rate = 0.02; % Baseline repair probability
    repair_strength = 0.03; % Weight recovery factor
    
    %% Initialize Storage
    time_points = 0:dt:t_max;
    n_steps = length(time_points);
    
    % Metrics storage for all realizations
    amyloid = zeros(n_steps, n_realizations);
    lscc = zeros(n_steps, n_realizations);
    efficiency = zeros(n_steps, n_realizations);
    repair_events = zeros(n_steps, n_realizations);
    V = zeros(n_steps, 1);
    %% Main Simulation
    for realization = 1:n_realizations
        fprintf('Realization %d/%d...\n', realization, n_realizations);
        
        % Create small-world network
        G = wattsStrogatz(N, k, p_rewire);
        G = make_directed(G, W_0);
        
        % Initialize variables
        A = 4; % Initial amyloid level
        
        for i = 1:n_steps
            %% Amyloid Dynamics (with microglial clearance)
            A = A + dt * (k_A * A * (1 - A / A_max));
            A = A - dt * clearance_rate * A; %*(A > A_th);
            A = max(0, A); % Prevent negative values
            amyloid(i, realization) = A; % Store for each realization
            
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
            repair_events(i, realization) = n_repairs; % Store for each realization
            
            %% Compute Network Metrics
            [~, sizes] = conncomp(G);
            lscc(i, realization) = max(sizes) / N; % Store for each realization
            
            % Global efficiency
            adj_mat = full(adjacency(G, 'weighted'));
            efficiency(i, realization) = global_efficiency_fun(adj_mat); % Store for each realization
        end
    end
    %% Calculate V for each N
        mean_lscc = mean(lscc(:,:), 2);
        mean_sc = mean(lscc(:,:).^2, 2);
        V = (mean_sc ./ mean_lscc) - mean_lscc;
        
        % Find critical time
        [~, tc_idx] = max(V);
        TC = time_points(tc_idx)
%         fprintf('Critical time T_c', TC);
    
    %% Calculate Means and Standard Deviations
    mean_amyloid = mean(amyloid, 2);
    std_amyloid = std(amyloid, 0, 2);
    mean_lscc = mean(lscc, 2);
    std_lscc = std(lscc, 0, 2);
    mean_efficiency = mean(efficiency, 2);
    std_efficiency = std(efficiency, 0, 2);
    mean_repairs = mean(repair_events, 2);
    std_repairs = std(repair_events, 0, 2);
    lscc_smooth = smoothdata(mean_lscc, 'sgolay', 21);
    d2 = gradient(gradient(lscc_smooth));
    [~, tc_idx] = max(d2);
    tc_in = time_points(tc_idx);
   disp(tc_in);
   figure,
   plot(time_points,V)
    %% Plot Results with Error Bars
    plot_metrics(time_points, mean_amyloid, std_amyloid, mean_lscc, std_lscc, mean_efficiency, std_efficiency, mean_repairs, std_repairs);
    
end

%% Plotting Function
function plot_metrics(time, mean_amyloid, std_amyloid, mean_lscc, std_lscc, mean_efficiency, std_efficiency, mean_repairs, std_repairs)
    % Create figure with 4 subplots
    figure('Position', [100, 100, 1000, 800]);
    
    % Amyloid levels
    subplot(4,1,1);
    errorbar(time, mean_amyloid, std_amyloid, 'b-', 'LineWidth', 2);
    ylabel('Amyloid Level');
    title('Amyloid Dynamics');
    grid on;
    
    % LSCC fraction
    subplot(4,1,2);
    errorbar(time, mean_lscc, std_lscc, 'r-', 'LineWidth', 2);
    ylabel('LSCC Fraction');
    title('Network Connectivity');
    grid on;
    
    % Global efficiency
    subplot(4,1,3);
    errorbar(time, mean_efficiency, std_efficiency, 'g-', 'LineWidth', 2);
    ylabel('Global Efficiency');
    title('Network Efficiency');
    grid on;
    
    % Repair events
    subplot(4,1,4);
    errorbar(time, mean_repairs, std_repairs, 'm-', 'LineWidth', 2);
    xlabel('Time (months)');
    ylabel('Repair Events');
    title('Network Repair Activity');
    grid on;
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