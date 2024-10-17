clc
clear
close all

WITH_CBFS = false;

% simulation parameters
DT = 0.01;
T = 1;
l_f = 1.4;
l_r = 1.6;

[f, g, full_dyn] = get_kinematic_bicycle_rearaxle(l_f, l_r);

% initialize system state
x0 = [-8; -8; 0; 2; 0];
t0 = 0;
x = x0;
% initial state of chain of integrators (dynamics of flat outputs)
x_df = [x0(1:2);
    x0(4) * cos(x0(3));
    x0(4) * sin(x0(3));
    zeros(2,1)];
x1 = x_df(1:2);
x2 = x_df(3:4);
x3 = x_df(5:6);

% tracking controller
rf = [8; 6];
r0 = x1 + (rf - x1) / 10;
r = r0;
phi = @(x1, x2, x3, r) 2 / T ^ 2 * (r - x1) - 2 / T * x2 - x3;
ALPHA = 5;

% safety controller
GAMMA = 1;
GAMMA_P = 1;
GAMMA_PP = 1;
o = [2; 2];
D = 2;
[h, h_dot, h_dot_dot, h_p, h_p_dot, h_pp, zeta] = build_parametric_functions(o, D, GAMMA, GAMMA_P, GAMMA_PP, ALPHA, phi);
if h(x1, x2, x3) < 0
    error('unsafe start')
end
if h_dot(x1, x2, x3) < 0
    GAMMA = -1 * h_dot(x1, x2, x3) / h(x1, x2, x3) * 1e1;
else
    GAMMA = 0.1;
end
if h_p_dot(x1, x2, x3) < 0
    GAMMA_P = -1 * h_p_dot(x1, x2, x3) / h_p(x1, x2, x3) * 1e1;
else
    GAMMA_P = 0.1;
end
GAMMA_PP = 5e-1;
[h, h_dot, h_dot_dot, h_p, h_p_dot, h_pp, zeta] = build_parametric_functions(o, D, GAMMA, GAMMA_P, GAMMA_PP, ALPHA, phi);

% log
logger.r = [];
logger.x = [];
logger.u = [];

% plot
h_br = plot_bicycle_rearaxle(x, l_f, l_r, []);
handle_x = scatter(x1' * [1; 0], x1' * [0; 1], 1000, 'r.');
handle_x_traj = plot(x1' * [1; 0], x1' * [0; 1], 'r', 'LineWidth', 2);
scatter(rf' * [1; 0], rf' * [0; 1], 1000, 'g.');
handle_r = scatter(nan, nan, 100, 'bo');
handle_r_traj = plot(nan, nan, 'b', 'LineWidth', 2);
handle_o = scatter(o' * [1; 0], o' * [0; 1], 1000, 'm.');
handle_o_radius = plot(o' * [1; 0] + D * cos(0 : 0.1 : 2*pi+0.1), ...
    o' * [0; 1] + D * sin(0 : 0.1 : 2*pi+0.1), 'm', 'LineWidth', 2);

% main simulation loop
for t = t0 : DT : norm(rf - r0)

    % compute reference
    r_dot = 1 * (rf - r) / norm(rf - r);
    r = r + r_dot * DT;

    % get states
    x_df(1:4) = [x(1:2);
        x(4) * cos(x(3));
        x(4) * sin(x(3))];

    x1 = x_df(1:2);
    x2 = x_df(3:4);
    x3 = x_df(5:6);

    % compute safety controller (closed form)
    zeta_t = zeta(x1, x2, x3, r);
    if zeta_t >= 0
        w = zeros(2,1);
    else
        w = -zeta_t * (x1 - o) / norm(x1 - o)^2;
    end

    % simulate flat output dynamics
    u_df = ALPHA * phi(x1, x2, x3, r) + ...
        WITH_CBFS * w;
    x_dot_df = [zeros(2), eye(2), zeros(2);
        zeros(2), zeros(2), eye(2);
        zeros(2), zeros(2), zeros(2)] * x_df + [zeros(2); zeros(2); eye(2)] * u_df;
    x_df = x_df + x_dot_df * DT;

    % compute system states and inputs from flat outputs
    x1 = x_df(1:2);
    x2 = x_df(3:4);
    x3 = x_df(5:6);

    if norm(x2) < eps
        a = 0;
        delta_dot = 0;
    else
        a = x2' * x3 / norm(x2)^2;
        delta_dot = (l_f + l_r) / norm(x2)^6 * ( ...
            (u_df(2) * x2(1) - u_df(1) * x2(2))*norm(x2)^3 - ...
            (x3(2) * x2(1) - x3(1) * x2(2))*3*norm(x2)*(x2'*x3) ...
            );
    end

    % simulate system dynamics
    u = [a; delta_dot];
    x_dot = f(x) + g(x) * u;
    x = x + x_dot * DT;

    % log
    logger.r(:, end+1) = r;
    logger.x(:, end+1) = x;
    logger.u(:, end+1) = u;

    % plot
    h_br = plot_bicycle_rearaxle(x, l_f, l_r, h_br);
    handle_x.XData = x(1);
    handle_x.YData = x(2);
    handle_x_traj.XData = [handle_x_traj.XData, x(1)];
    handle_x_traj.YData = [handle_x_traj.YData, x(2)];
    handle_r.XData = r(1);
    handle_r.YData = r(2);
    handle_r_traj.XData = [handle_r_traj.XData, r(1)];
    handle_r_traj.YData = [handle_r_traj.YData, r(2)];
    drawnow limitrate

end