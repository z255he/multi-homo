function dydt = ode3(t, y, param)
    R0 = 10;
    %R0 = param.R0;
    gamma = 36;
    mu = 1/60;
    beta = R0 * (gamma + mu);
    epsilon = 0.9;
    kappa = 1000;
    delta = 8e-4;
    lambda = [0; 0; 0];
    N = 3;
    
    s = y(1);
    i = y(2);
    v = y(3);
    xi = y(4:6);
    x = param.pop * xi / sum(param.pop);
    
    omega = 3e-4;
    
    %EP = -omega + delta * (param.h' .* xi + (1 - param.h') * sum(xi));
    %EA = -i + delta * (param.h' .* (1 - xi) + (1 - param.h') * (N - sum(xi)));
    EP = -omega + delta * (xi + (1 - param.h')./(sum(param.pop) - param.pop') .* (param.pop * xi - param.pop' .* xi));
    EA = -i + delta * ((1 - xi) + (1 - param.h')./(sum(param.pop) - param.pop') .* (param.pop * (1 - xi) - param.pop' .* (1 - xi)));

    EAtoP = max(EP' - EA, 0);
    EPtoA = max(EA' - EP, 0);
    K = kappa*param.pop/sum(param.pop);
    
    M = (K .* param.H) .* ((eye(N) - diag(xi)) * EAtoP * diag(xi) - diag(xi) * EPtoA * (eye(N) - diag(xi))) + lambda .* (1 - xi);
    
    f_s = mu * (1 - epsilon * x) - beta * s * i - mu * s;
    f_i = beta * s * i - (gamma + mu) * i;
    f_v = mu * epsilon * x - mu * v;
    f_xi = M * ones(N, 1);
    
    dydt = [f_s; f_i; f_v; f_xi];
end