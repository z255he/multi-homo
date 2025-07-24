function dydt = ode3tv(t, y, param)
    R0 = 10;
    gamma = 24;
    mu = 1/60;
    beta = R0 * (gamma + mu);
    epsilon = 0.9;
    kappa = 1000;
    delta = 8e-4;
    N = 3;
    
    s = y(1);
    i = y(2);
    v = y(3);
    xi = y(4:6);
    x = param.pop * xi / sum(param.pop);
    
    omega = 3e-4;
    
    % Logistic decrease
    %L = 0.1;
    %ht = 2 * (param.h - L)/(1 + exp(param.alpha * t)) + L;

    % Logistic increase
    L = 0.9;
    ht = 2 * (L - param.h)/(1 + exp(-param.alpha * t)) + 2 * param.h - L;

    % Linear decrease
    %L = 0.1;
    %ht = max(-param.alpha * t + param.h, L);

    % Linear increase
    %L = 0.9;
    %ht = min(param.alpha * t + param.h, L);

    % Piecewise constant disturbance
    %t_start = 30;
    %t_end = 60;
    %H_window = double(t >= t_start) .* double(t <= t_end);
    %A = 0.1;
    %ht = max(min(param.h + A * H_window, 0.9), 0.1);

    % Pseudo-random perturbation
    %perturb = sum(param.C .* sin(param.w * t + param.phi));
    %ht = max(min(param.h + perturb, 0.9), 0.1);

    % Homophily matrix
    H = zeros(3); % Initialize n×n matrix

    for j = 1:3
        H(j, :) = (1 - ht(j)); % Fill row i with (1 - h_i)
        H(j, j) = 1 + 2 * ht(j); % Set diagonal element
    end
    H = 1/3 * H;
    
    %EP = -omega + delta * (ht' .* xi + (1 - ht') * sum(xi));
    %EA = -i + delta * (ht' .* (1 - xi) + (1 - ht') * (N - sum(xi)));
    EP = -omega + delta * (xi + (1 - param.h')./(sum(param.pop) - param.pop') .* (param.pop * xi - param.pop' .* xi));
    EA = -i + delta * ((1 - xi) + (1 - param.h')./(sum(param.pop) - param.pop') .* (param.pop * (1 - xi) - param.pop' .* (1 - xi)));

    EAtoP = max(EP' - EA, 0);
    EPtoA = max(EA' - EP, 0);
    K = kappa*param.pop/sum(param.pop);
    
    M = (K .* H) .* ((eye(N) - diag(xi)) * EAtoP * diag(xi) - diag(xi) * EPtoA * (eye(N) - diag(xi)));
    
    f_s = mu * (1 - epsilon * x) - beta * s * i - mu * s;
    f_i = beta * s * i - (gamma + mu) * i;
    f_v = mu * epsilon * x - mu * v;
    f_xi = M * ones(N, 1);
    
    dydt = [f_s; f_i; f_v; f_xi];
end