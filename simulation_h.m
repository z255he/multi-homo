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
tf = 500;  % 2000 for infty
t_uniform = linspace(0, tf, 10000)';

% Initial conditions
y0 = [0.01; 1e-6; 0.8; 0.69; 0.95; 0.2];

% Fixed parameters
param_base.h = [0.5, 0.8, 0.8];
param_base.pop = [8e4, 1e4, 1e4];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% Parameter range for h
h_values = 0.1:0.01:0.9;
%r0_values = [0.8, 4, 8, 10];
%num_sims = length(h_values);
num_sims = length(h_values)*length(h_values);
%num_sims = length(h_values)*length(r0_values);

% Preallocate result storage (cell array in case y1 has multiple columns)
results = cell(num_sims, 1);

%{
for i = 1:num_sims
    % Update parameter
    param = param_base;
    param.h(1) = h_values(i);
    param.H = homophily_matrix(param.h);  % Assuming this returns a 3x3 matrix

    % Solve the ODE
    [t_sol, y_sol] = ode45(@(t, y) ode3(t, y, param), [0, tf], y0, options);

    % Interpolate to uniform time grid
    y_interp = interp1(t_sol, y_sol, t_uniform, 'linear');

    % Store result
    results{i} = struct('h', h_values(i), 't', t_uniform, 'y', y_interp);
end

% Define a colormap
cmap = [...
    240, 249, 232;
    204, 235, 197;
    168, 221, 181;
    123, 204, 196;
     78, 179, 211;
     43, 140, 190;
      8,  88, 158] / 255;
cmap_interp = interp1(linspace(0,1,size(cmap,1)), cmap, linspace(0,1,num_sims));

figure; hold on;
for i = 1:num_sims
    t = results{i}.t;
    y = results{i}.y(:, 4);
    h = h_values(i);

    x_vals = h * ones(size(t));
    y_vals = t;
    z_vals = y;

    plot3(x_vals, y_vals, z_vals, 'Color', cmap_interp(i, :), 'LineWidth', 1.2);
end

xlabel('Homophily');
ylabel('Time');
zlabel('Pro-vaccine');
view(60, 30);
grid on;

figure; hold on;
for i = 1:num_sims
    t = results{i}.t;
    y = results{i}.y(:, 6);
    h = h_values(i);

    x_vals = h * ones(size(t));
    y_vals = t;
    z_vals = y;

    plot3(x_vals, y_vals, z_vals, 'Color', cmap_interp(i, :), 'LineWidth', 1.2);
end

xlabel('Homophily');
ylabel('Time');
zlabel('Infected');
view(60, 30);
grid on;

% Extract number of time points and parameter values
num_time_points = length(t_uniform);

% Preallocate matrix to store interpolated values of y(:,1)
Z = zeros(num_time_points, num_sims);  % rows: time, cols: parameter

for i = 1:num_sims
    Z(:, i) = results{i}.y(:, 4);  % extract the first variable of the solution
end

% Create meshgrid for plotting
[H_grid, T_grid] = meshgrid(h_values, t_uniform);  % X: parameter, Y: time

% Define a colormap
cmap = [...
     13,   8, 135;
     75,  3, 161;
    125,  3, 168;
    168, 34, 150;
    203, 70, 121;
    234,121,  83;
    252,185,  44;
    240,249,  33] / 255;

% Plot
figure;
surf(T_grid, H_grid, Z, 'EdgeColor', 'none');
xlabel('Time');
ylabel('Homophily');
zlabel('Pro-vaccine');
colormap(interp1(linspace(0,1,size(cmap,1)), cmap, linspace(0,1,100)));
colorbar;
view(0, 90);  % Adjust viewing angle
%}

%{
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
        labels{i} = sprintf('h_%.2f', results{i}.h);
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
%}

%
% Vary both h2 and h3
for i = 1:length(h_values)
    for j = 1:length(h_values)
        % Update parameter
        param = param_base;
        %param.R0 = r0_values(i);
        param.h(2) = h_values(i);
        param.h(3) = h_values(j);
        param.H = homophily_matrix(param.h);  % Assuming this returns a 3x3 matrix

        % Solve the ODE
        [t_sol, y_sol] = ode45(@(t, y) ode3(t, y, param), [0, tf], y0, options);

        % Interpolate to uniform time grid
        y_interp = interp1(t_sol, y_sol, t_uniform, 'linear');

        % Store result
        idx = (i - 1)*length(h_values) + j;
        results{idx} = struct('v1', h_values(i), 'v2', h_values(j), 't', t_uniform, 'y', y_interp);
    end
end

% Save files
num_vars = size(results{1}.y, 2);         % number of variables in y_interp
num_time = length(results{1}.t);          % number of time points
num_sims = length(results);               % number of parameter values

for var_idx = 1:num_vars
    % Initialize matrix: first column is time
    data_matrix = zeros(num_time, num_sims + 1);
    data_matrix(:, 1) = results{1}.t;     % assume same time vector for all
    labels = cell(1, num_sims);

    % Fill in solution data
    for i = 1:num_sims
        data_matrix(:, i+1) = results{i}.y(:, var_idx);
        labels{i} = sprintf('h2_%.2f_h3_%.2f', results{i}.v1, results{i}.v2);
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
%}