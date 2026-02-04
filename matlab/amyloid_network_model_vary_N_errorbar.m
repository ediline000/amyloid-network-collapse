function amyloid_network_model_vary_N_errorbar
    clc
    %% Network sizes to test
    N_values = [8000,10000];

    %% Other parameters
    k = 20;
    p_rewire = 0.1;
    t_max = 500;
    dt = 1;
    n_realizations = 20;

    A_max = 150;
    k_A = 0.02;
    A_th = 25;
    clearance_rate = 0.005;

    W_0 = 1;
    W_min = 0.5*W_0 * exp(-1);
    tau_synapse = 120;

    base_repair_rate = 0.02;
    repair_strength = 0.03;

    time_points = 0:dt:t_max;
    n_steps = length(time_points);

    %% Initialize storage for critical times & LSCC curves
    critical_times_d2_all = zeros(length(N_values), n_realizations);
    critical_times_2nd_all = zeros(length(N_values), n_realizations);

    mean_lscc_allN = cell(length(N_values), 1);
    std_lscc_allN = cell(length(N_values), 1);
    mean_second_lscc_allN = cell(length(N_values), 1);
    std_second_lscc_allN = cell(length(N_values), 1);

    %% Loop over network sizes
    for idx = 1:length(N_values)
        N = N_values(idx);
        fprintf('\n=== Running for N = %d ===\n', N);

        lscc_all = zeros(n_realizations, n_steps);
        second_lscc_all = zeros(n_realizations, n_steps);

        for realization = 1:n_realizations
            fprintf('Realization %d/%d...\n', realization, n_realizations);

            lscc = zeros(1, n_steps);
            second_lscc = zeros(1, n_steps);
            A = 4;

            G = wattsStrogatz(N, k, p_rewire);
            G = make_directed(G, W_0);

            for i = 1:n_steps
                %% Amyloid dynamics
                A = A + dt*(k_A*A*(1 - A/A_max));
                A = A - dt*clearance_rate*A;
                A = max(0, A);

                %% Synaptic degradation
                p_toxicity = 1 - exp(-0.001*max(0, A - A_th));
                edges_to_remove = [];
                for e = 1:numedges(G)
                    if rand() < p_toxicity
                        G.Edges.Weight(e) = W_0 * exp(-(i*dt)/tau_synapse);
                    end
                    if G.Edges.Weight(e) < W_min
                        edges_to_remove = [edges_to_remove; e];
                    end
                end
                G = rmedge(G, edges_to_remove);

                %% Repair
                current_repair_rate = base_repair_rate * (1 - exp(-0.01*(A_th - A))) / sqrt(N);
                [G, ~] = repair_network(G, W_0, current_repair_rate, repair_strength);

                %% Compute LSCC and second LSCC
                [~, sizes] = conncomp(G);
                sorted_sizes = sort(sizes, 'descend');
                lscc(i) = sorted_sizes(1) / N;
                if numel(sorted_sizes) >= 2
                    second_lscc(i) = sorted_sizes(2) / N;
                else
                    second_lscc(i) = 0;
                end
            end

            lscc_all(realization, :) = lscc;
            second_lscc_all(realization, :) = second_lscc;

            %% Compute critical times
            lscc_smooth = smoothdata(lscc, 'sgolay', 21);
            d2 = gradient(gradient(lscc_smooth));
            [~, tc_idx_d2] = max(d2);
            critical_times_d2_all(idx, realization) = time_points(tc_idx_d2);

            [~, tc_idx_2nd] = max(second_lscc);
            critical_times_2nd_all(idx, realization) = time_points(tc_idx_2nd);
        end

        %% Compute mean & std over realizations
        mean_lscc = mean(lscc_all, 1);
        std_lscc = std(lscc_all, [], 1);

        mean_second_lscc = mean(second_lscc_all, 1);
        std_second_lscc = std(second_lscc_all, [], 1);

        %% Store for later plotting
        mean_lscc_allN{idx} = mean_lscc;
        std_lscc_allN{idx} = std_lscc;
        mean_second_lscc_allN{idx} = mean_second_lscc;
        std_second_lscc_allN{idx} = std_second_lscc;
    end

    %% ✅ Plot LSCC and 2nd LSCC over time with error bars for all N in same figure
    figure;
    colors = lines(length(N_values));

    subplot(2,1,1); hold on;
    for idx=1:length(N_values)
        errorbar(time_points, mean_lscc_allN{idx}, std_lscc_allN{idx}, ...
            'Color', colors(idx,:), 'LineWidth', 1.5, 'DisplayName', ['N=' num2str(N_values(idx))]);
    end
    xlabel('Time (months)'); ylabel('LSCC fraction');
    title('LSCC over time with error bars'); legend show; grid on;

    subplot(2,1,2); hold on;
    for idx=1:length(N_values)
        plot(time_points, mean_second_lscc_allN{idx}, ...
            'Color', colors(idx,:), 'LineWidth', 1.5, 'DisplayName', ['N=' num2str(N_values(idx))]);
    end
    xlabel('Time (months)'); ylabel('2nd LSCC fraction');
    title('2nd LSCC over time with error bars'); legend show; grid on;

    %% Compute mean & std of critical times over realizations
    mean_ct_d2 = mean(critical_times_d2_all, 2);
    std_ct_d2 = std(critical_times_d2_all, [], 2);
    mean_ct_2nd = mean(critical_times_2nd_all, 2);
    std_ct_2nd = std(critical_times_2nd_all, [], 2);

    %% Plot critical times vs N with error bars
    figure;
    plot(N_values, mean_ct_d2, 'bo-', 'LineWidth', 1.5, 'DisplayName','Curvature');
    hold on;
    plot(N_values, mean_ct_2nd, 'ms-', 'LineWidth', 1.5, 'DisplayName','2nd LSCC peak');
    xlabel('Network size N'); ylabel('Critical time (months)');
    title('Critical time vs Network size'); legend show; grid on;

    %% Save results
    T = table(N_values(:), mean_ct_d2, std_ct_d2, mean_ct_2nd, std_ct_2nd, ...
        'VariableNames', {'N', 'MeanCt_D2', 'StdCt_D2', 'MeanCt_SecondLSCC', 'StdCt_SecondLSCC'});
    writetable(T, 'critical_times_results2.csv');
    save('critical_times_results2.mat', 'N_values', 'critical_times_d2_all', 'critical_times_2nd_all', ...
        'mean_ct_d2', 'std_ct_d2', 'mean_ct_2nd', 'std_ct_2nd');

    fprintf('\n✅ Results saved to critical_times_results2.mat and critical_times_results2.csv\n');
end

%% Helper functions (same as yours, copy as-is)
function G = wattsStrogatz(N, k, beta)
    G = graph(); G = addnode(G, N);
    for i = 1:N
        for j = 1:k/2
            neighbor = mod(i + j - 1, N) + 1;
            G = addedge(G, i, neighbor, 1);
        end
    end
    edges = table2array(G.Edges);
    for e = 1:size(edges,1)
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
    n_repairs = 0;
    for e=1:numedges(G)
        if rand() < p_repair
            G.Edges.Weight(e) = min(W_0, G.Edges.Weight(e)*(1+repair_strength));
            n_repairs = n_repairs+1;
        end
    end
    if rand() < p_repair/2
        nodes = randperm(numnodes(G),2);
        if ~findedge(G, nodes(1), nodes(2))
            G = addedge(G, nodes(1), nodes(2), W_0*0.5);
            n_repairs = n_repairs+1;
        end
    end
end
