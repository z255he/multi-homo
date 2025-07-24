function H = homophily_matrix(h)
    n = length(h);
    H = zeros(n); % Initialize n×n matrix

    for i = 1:n
        H(i, :) = (1 - h(i)); % Fill row i with (1 - h_i)
        H(i, i) = 1 + (n - 1) * h(i); % Set diagonal element
    end
    H = 1/n * H;
end

% Final time and uniform time vector
tf = 200;
t_uniform = linspace(0, tf, 10000)';

% Initial conditions
y0 = [0.01; 1e-6; 0.8; 0.65; 0.95; 0.2];

% Parameters for standard case
param.h = [0.5, 0.8, 0.8];
param.pop = [8e4, 1e4, 1e4];
param.H = homophily_matrix(param.h);

% Monotonic logistic trend
alpha = [0, 0.02, 0.05];
num_sims = length(alpha);

% Simulation
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
for i = 1:num_sims
    param.alpha = alpha(i);
    [t2, y2] = ode45(@(t, y) ode3tv(t, y, param), [0, tf], y0, options);
    y_interp = interp1(t2, y2, t_uniform, 'linear');
    results{i} = struct('alpha', alpha(i), 't', t_uniform, 'y', y_interp);
end

% Save files
num_vars = size(results{1}.y, 2);         % number of variables in y_interp
num_time = length(results{1}.t);          % number of time points
num_sims = length(results);               % number of parameter values

for var_idx = 1:num_vars
    % Initialize matrix: first column is time
    data_matrix = zeros(num_time, num_sims + 1);
    data_matrix(:, 1) = results{1}.t;     % assume same time vector for all

    % Fill in solution data
    for i = 1:num_sims
        data_matrix(:, i+1) = results{i}.y(:, var_idx);
        labels{i} = sprintf('alpha_%.2f', results{i}.alpha);
    end

    % Combine with 'time' as the first header
    header = [{'time'}, labels];

    % Write to CSV
    filename = sprintf('solution_var%d.csv', var_idx);
    fid = fopen(filename, 'w');
    fprintf(fid, '%s,', header{1:end-1});
    fprintf(fid, '%s\n', header{end});
    fclose(fid);
    dlmwrite(filename, data_matrix, '-append');
end